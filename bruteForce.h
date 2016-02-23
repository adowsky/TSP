#ifndef TSP_BRUTEFORCE_H
#define TSP_BRUTEFORCE_H

#include "Algorithms.h"
#include "Graph.h"

class bruteForce : public Algorithms {
public:
    bruteForce(Graph *graph);
    void runAlgorithm();
    void printResult();

private:
    Graph* graph;

    void initialize();
    int* createSortedSequence();
    void checkAllResults(int* sortedSequence);
    void setBestCost(int* sequence, int cost);
};


#endif //TSP_BRUTEFORCE_H
