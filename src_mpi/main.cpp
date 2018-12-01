#include <iostream>
#include <mpi.h>

#include "map.h"
#include "debug.h"
#include "dijkstra.h"

//note; all LOG(x) calls are only defined if compiled with -DDEBUG flax in makefile

int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cout << "Usage: " << argv[0] << " <testcase file>" << std::endl;
        return -1;
    }
    Map m = Map::fromFile(argv[1]);
    m.printWeights();

    int mpiNodesCount, mpiNodeId;
    const int mpiRootId = 0;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiNodesCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiNodeId);
    std::cout << "Hello from " << mpiNodeId << std::endl;

    if (mpiNodeId == mpiRootId) {
        Map m = Map::fromFile(argv[1]);
        
        auto n = m.getNodesNames();
        auto initialNodeName = *n.begin();
        auto goalNodeName = *(n.end()-1);

        //a lot of this should be optimized; printWeights is a O(n^2) function which adds a big constant factor...
        
        std::cout << "Dijkstra search algorithm" << std::endl;
        std::cout << "Starting at node: " << initialNodeName << std::endl;
        std::cout << "Ending at node: " << goalNodeName << std::endl;
        std::cout << "Consecutive nodes (A, B, ...) weights: " << std::endl;
        m.printWeights();

        std::cout << "Searching..." << std::endl;
        dijkstra(m, initialNodeName, goalNodeName, mpiNodesCount);
    }
    else
        dijkstraWorker(mpiNodeId, mpiNodesCount);

    MPI_Finalize();

    return 0;
}
