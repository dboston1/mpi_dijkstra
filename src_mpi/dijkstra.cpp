#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <limits>

#include <mpi.h>

#include "debug.h"
#include "dijkstra.h"
#include "map.h"

//added for std::sqrt() functionality
#include <math.h>

static const auto INF = std::numeric_limits<int>::max();
const int mpiRootId = 0;


//this function takes nodesCount, mpiNodesCount (num of processors - host), and mpiNodeId
// it essentially figures out which processor calculated the distance from "fromNode" to "toNode"
std::pair<int, int> getMpiWorkerNodeRanges(int nodesCount, int mpiNodesCount, int mpiNodeId) {
    mpiNodesCount -= 1; // counting only worker nodes
    auto fromNode = (nodesCount / mpiNodesCount) * (mpiNodeId - 1);
    auto toNode = (nodesCount / mpiNodesCount) * mpiNodeId - 1;

    auto restNodes = nodesCount % mpiNodesCount;

    if (mpiNodeId - 1 < restNodes) {
        fromNode += mpiNodeId - 1;
        toNode += mpiNodeId;
    }
    else {
        fromNode += restNodes;
        toNode += restNodes;
    }

    return std::pair<int, int>(fromNode, toNode);
}
    
    
auto isNeighbour(auto currentNode, auto node, auto dim){
    int sourceNode = dim*dim+1;
    if(currentNode == 0 && node == sourceNode)
       return false;
    if(currentNode == 0){
        return ((node-1) % dim) == 0;
    }
    if(node == sourceNode){
       return (currentNode % dim) == 0;
    }

    int i_curr = (int)((currentNode-1) / dim);
    int j_curr = (currentNode-1) % dim;
    int i_to = (int)((node-1) / dim);
    int j_to = (node-1) % dim;

    if(i_curr == 0){
        if((i_curr == i_to) || (i_curr == (i_to - 1))){
            return(j_curr == (j_to -1));
        }
    }
    else if(i_curr == (dim -1)){
        if((i_curr == i_to) || (i_curr == (i_to + 1))){
            return(j_curr == (j_to -1));
        }
    }
    else{
        if(j_curr == (j_to - 1)){
            if((i_curr == i_to) || (i_curr == (i_to - 1)) || (i_curr == (i_to + 1))){
                return true;
            }
        }
    }
    return false;
}  
    

//only the "main" processer (id == 0) uses this function
void dijkstra(const Map& m, const std::string& initialNodeName, const std::string& goalNodeName, const int mpiNodesCount) {
    
    const auto& weights = m.getWeights();
    
    const auto& nodesNames = m.getNodesNames();
    auto nodesCount = nodesNames.size();
    auto dim = std::sqrt(nodesCount-2);

    std::vector<int> distances(nodesCount);
    std::vector<int> prevNodes(nodesCount);

    std::set<int> visited;
    
    // 0u is an unsigned short; has range [0, 65k] so should be fine for nodesCount <= 65k
    for(auto node=0u; node<nodesCount; ++node) {
        distances[node] = INF;
        prevNodes[node] = INF;
    }

    auto indexOf = [&] (auto nodeName) { return std::stoi(nodeName)-1; };
        
    auto isVisited = [&] (auto node) { return visited.find(node) != visited.end(); };

    auto initialNode = static_cast<int>(0);
    auto currentNode = initialNode;
    auto goalNode = static_cast<int>(std::stoi(goalNodeName));

    int workerNodes = mpiNodesCount - 1;
    int nodesStep = nodesCount / workerNodes;

    
    // so here, we broadcast twice;
    // first, the nodesCount, initialNode name, and goalNode name (all as ints)
    // second, the edge weights from the initialNode to every other node (i.e. the first column of the adjacency matrix, since [i][j] is non-negative if 
    //          the edge (i,j) exists 
    
    std::cout << "Sending initial data to workers..." << std::endl;
    int data[3] = {nodesCount, initialNode, goalNode};
    MPI_Bcast(&data, 3, MPI_INT, mpiRootId, MPI_COMM_WORLD);

    //This is sending a pointer to each row in the adjacency matrix...
    for(auto i=0u; i<dim; ++i)
        MPI_Bcast((int*)&weights[i][0], dim, MPI_INT, mpiRootId, MPI_COMM_WORLD);
    distances[initialNode] = 0;
    
    
    
    while (1) {
        std::cout << "Sending currNode=" << currentNode << " distance=" << distances[currentNode] << std::endl;
        // send current node (starting with initialNode) and its value in distances to all other processors
        int data[2] = {currentNode, distances[currentNode]};
        MPI_Bcast(&data, 2, MPI_INT, mpiRootId, MPI_COMM_WORLD);

        // receive distances calculated by workers
        // for each processor (besides host):
        for(auto mpiNodeId=1; mpiNodeId<mpiNodesCount; ++mpiNodeId) {
            // see function at top of file for getMPIWorkerNodeRanges()
            const auto nodeRanges = getMpiWorkerNodeRanges(nodesCount, mpiNodesCount, mpiNodeId);
            const auto fromNode = nodeRanges.first;
            const auto toNode = nodeRanges.second;
            MPI_Recv(&distances[fromNode], toNode - fromNode + 1, MPI_INT, mpiNodeId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            std::cout << "Recv from " << mpiNodeId << std::endl;
        }

        // test for goal
        if (currentNode == goalNode) {
            std::cout << "Goal node found" << std::endl;
            for(auto mpiNodeId=1; mpiNodeId<mpiNodesCount; ++mpiNodeId) {
                // see function at top of file for getMPIWorkerNodeRanges()
                const auto nodeRanges = getMpiWorkerNodeRanges(nodesCount, mpiNodesCount, mpiNodeId);
                const auto fromNode = nodeRanges.first;
                const auto toNode = nodeRanges.second;
                MPI_Recv(&prevNodes[fromNode], toNode - fromNode + 1, MPI_INT, mpiNodeId, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            
            LOGv(distances);
            LOGv(prevNodes);

            std::vector<int> stack;
            while(currentNode != initialNode) {
                stack.push_back(currentNode);

                auto prev = prevNodes[currentNode];
                currentNode = prev;
            }

            LOGv(stack);
            std::cout << "Total cost: " << distances[goalNode] << std::endl;
            std::cout << "Path: " << nodesNames[currentNode];

            for(auto it=stack.rbegin(); it != stack.rend(); ++it) {
                auto nextNodeName = nodesNames[*it];
                std::cout << " -> " << nextNodeName;
            }
            std::cout << std::endl;

            break;
        }
        
        //otherwise, goal has not been found yet:
        LOG("Visited " << currentNode);
        // mark current node as visited
        visited.insert(currentNode);

        auto minCost = std::numeric_limits<int>::max(); 
        auto nextNode = -1;
        
        
        //this part is the "greedy" part of algorithm: find next unvisited node with minimal path from source to it and set it 
        // as next node. If none exists, it will remain -1 at end, at which point program will exit after printing "Path not Found"
        for(auto node=0u; node<nodesCount; ++node) {
            auto totalCost = distances[node];

            if (!isVisited(node) && totalCost < minCost) {
                minCost = totalCost;
                nextNode = node;
            }
        }
        
        currentNode = nextNode;
        LOG("Next currentNode = " << currentNode);

        if (currentNode == -1) {
            std::cout << "Path not found" << std::endl;
            return;
        }
    }
}


// all processors except main use this function (only)
void dijkstraWorker(int mpiNodeId, int mpiNodesCount) {
    int data[3];
    // recieve nodesCount, intialNode, and goalNode from host processor
    MPI_Bcast(&data, 3, MPI_INT, 0, MPI_COMM_WORLD);

    int nodesCount = data[0];
    int dim = std::sqrt(nodesCount-2);
    int initialNode = data[1];
    int goalNode = data[2];

    std::cout << "mpiID=" << mpiNodeId << " bcast recv nodesCount=" << nodesCount << " initialNode=" << initialNode << " goalNode=" << goalNode << std::endl;

    std::vector<std::vector<int>> weights(nodesCount);
    std::vector<int> distances(nodesCount, INF);
    std::vector<int> prevNodes(nodesCount, INF);
    std::set<int> visited;

    auto isVisited = [&] (auto node) { return visited.find(node) != visited.end(); };
    
    //populate weights from host broadcast 
    for(auto i=0u; i<dim; ++i) {
        weights[i].resize(dim);
        MPI_Bcast((int*)&weights[i][0], dim, MPI_INT, mpiRootId, MPI_COMM_WORLD);
    }
    // nodeRanges is the nodes this processor will be checking; uses function at top of file
    const auto nodeRanges = getMpiWorkerNodeRanges(nodesCount, mpiNodesCount, mpiNodeId);
    const auto fromNode = nodeRanges.first;
    const auto toNode = nodeRanges.second;
    std::cout << "mpiID=" << mpiNodeId << " from=" << fromNode << " to=" << toNode << std::endl;

    // real work
    while (1) {
        std::cout << "~~~~" << std::endl;
        // get currentNode, and path distance from initialNode to currentNode (guaranteed to be shortest path seen so far)
        MPI_Bcast(&data, 2, MPI_INT, mpiRootId, MPI_COMM_WORLD);

        int currentNode = data[0];
        distances[currentNode] = data[1];
        std::cout << "mpiId=" << mpiNodeId << " bcast recv currNode=" << currentNode << " dist=" << distances[currentNode] << std::endl;

        // goal node not found
        if (currentNode == -1)
            return;

        for (auto node=fromNode; node<=toNode; ++node) {
            if (isVisited(node)) {
                LOG("Node " << node << " already visited - skipping.");
                continue;
            }

            if (isNeighbour(currentNode, node, dim)) {
                int row_index = (int)((node-1) / dim);
                int col_index = (node-1) % dim;
                auto nodeDistance = 0;
                if(node != (nodesCount -1)){
                    nodeDistance = weights[row_index][col_index];
                }
                
                
                auto totalCostToNode = distances[currentNode] + nodeDistance;
                LOG("Node " << node << " is neighbour of " << currentNode << " (distance: " << nodeDistance << ", totalCostToNode: " << totalCostToNode << ")");
                if (totalCostToNode < distances[node]) {
                    distances[node] = totalCostToNode;
                    prevNodes[node] = currentNode;
                    LOG("New total cost is less than the old, replacing");
                }
            }
        }

        visited.insert(currentNode);
        // send all calculated distances
        MPI_Send(&distances[fromNode], toNode - fromNode + 1, MPI_INT, mpiRootId, 0, MPI_COMM_WORLD);

        // goal node found
        if (currentNode == goalNode) {
            MPI_Send(&prevNodes[fromNode], toNode - fromNode + 1, MPI_INT, mpiRootId, 0, MPI_COMM_WORLD);
            return;
        }
        std::cout << "======" << std::endl;
    }
}
