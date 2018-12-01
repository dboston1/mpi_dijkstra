#pragma once

#include <string>

#include "map.h"

auto isNeighbour(auto currentNode, auto node);

void dijkstra(const Map& m, const std::string& initialNodeName, const std::string& goalNodeName, const int mpiNodesCount);

void dijkstraWorker(int mpiNodeId, int mpiWorkerNodes);
