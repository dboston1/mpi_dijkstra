#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>

#include "map.h"
#include "debug.h"

//added for sqrt of vertices count calculation
#include <math.h>

Map::Map(int dim) : weights(dim) {
    for(int i=0; i<dim; ++i) {
        weights[i].resize(dim);
    }
}


Map Map::fromFile(const std::string& str, const char delimiter) {
    return fromFile(std::ifstream(str.c_str()), delimiter);
}


Map Map::fromFile(std::ifstream&& istream, const char delimiter) {
    std::string header;

    if (std::getline(istream, header, ':') && header == "vertices" && std::getline(istream, header)) {
        int verticesCount = std::stoi(header);
        int dim = std::sqrt(verticesCount);
        Map m(dim);        

        for(auto i=0; i<dim; ++i) {
            std::string line;
            std::getline(istream, line); 

            std::stringstream linestream(line);
            for(auto j=0; j<dim; ++j) {
                std::string numstr;
                std::getline(linestream, numstr, delimiter);
                // 9 - the number because we want to find the LONGEST path, this conversion makes it a shortest paths problem 
                m.weights[i][j] = 9 - (std::stoi(numstr));
            }
        }
        std::stringstream ss;

        for(auto i=0; i<verticesCount+2; ++i) {
            ss << i;
            std::string nodeName = ss.str();
            ss.str(std::string());
            //std::cout << "Adding node: " << nodeName << std::endl;
            m.nodesNames.push_back(std::string(nodeName));
        }
        
        
        return m;
    }

    throw std::runtime_error("Wrong map file format");
}

void Map::printWeights() const {
    for(auto i=0u; i<getSize(); ++i) {
        for(auto j=0u; j<getSize(); ++j)
            std::cout << std::setw(4) << (9 - weights[i][j]) << " ";
        std::cout << std::endl;
    }
}
