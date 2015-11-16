#include <iostream>
#ifndef TSP_GRAPH_H
#define TSP_GRAPH_H
using namespace std;

class Graph {
public:
    Graph();
    Graph(string s);
    int** getGrpah();
    int getFirstVertex();
    int getVertices();
    void printGraph();
private:
    int** graph;
    int firstVertex;
    int vertices;

};

#endif //TSP_GRAPH_H
