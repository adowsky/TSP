#include "Graph.h"
#include <stdlib.h>
Graph::Graph() {
    cout<<"Vertices number: ";
    cin>>vertices;
    graph = new int*[vertices];
    srand(time(NULL));
    // creating arrays and filling them with random values.
    for (int i = 0; i < vertices; ++i) {
        graph[i] = new int[vertices];
    }
    for (int i = 0; i < vertices; ++i) {

        for (int j = 0; j < vertices; ++j) {
            if (j == i)
                graph[i][j] = 0;
            else if (j < i)
                graph[j][i] = graph[i][j];
            else
                graph[j][i] = rand() % (vertices * 10) + 1;
        }
    }
    firstVertex = rand()%vertices;
}
int** Graph::getGrpah() {return graph; }

int Graph::getFirstVertex() {return firstVertex; }

int Graph::getVertices() {return vertices; }
void Graph::printGraph() {
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            cout<<graph[i][j]<<" ";
        }
        cout<<endl;
    }
}