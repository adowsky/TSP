
#include "bruteForce.h"
#include <algorithm>

bruteForce::bruteForce(Graph *graph) {
    this->graph = graph;
    bestWay = NULL;
}

void bruteForce::initialize() {
    bestCost = 0;
    if(bestWay == NULL)
        bestWay = new int[graph->getVertices()];
    bestWay[0] = graph->getFirstVertex();
}
int* bruteForce::createSortedSequence() {
    int vertices = graph->getVertices();
    int* sequence= new int[vertices];
    for (int j = 0; j < vertices; ++j) {
        sequence[j] = j;
    }
    return sequence;
}

void bruteForce::setBestCost(int *sequence, int cost) {
    bestCost = cost;
    for (int i = 0; i < graph->getVertices(); ++i) {
        bestWay[i] = sequence[i];
    }
}
void bruteForce::checkAllResults(int *sortedSequence) {

    int worstCost = 0;
    int** wayGraph = graph->getGrpah();
    int vertices = graph->getVertices();
    int i=0;
    do{
        int cost = 0;
        for (int i = 0; i < vertices; ++i) {
            cost += wayGraph[sortedSequence[i]][sortedSequence[(i+1)%(vertices)]];

        }
        if(cost<bestCost||bestCost==0)
            setBestCost(sortedSequence, cost);
        if(cost>worstCost)
            worstCost= cost;
        i++;
    }while(next_permutation(sortedSequence, sortedSequence+vertices-1 ));
    cout<<"Worst Cost: "<<worstCost<<endl;
}
void bruteForce::runAlgorithm() {

    clock_t timer = clock();
    initialize();
    int* sequence = createSortedSequence();
    checkAllResults(sequence);
    algorithmTime = (clock()-timer)/(double)CLOCKS_PER_SEC;
}

void bruteForce::printResult() {
    int first = graph->getFirstVertex();
    int firstIndex = -1;
    bool found = false;
    int i=0;
    while(!found){
        if(bestWay[i] == first){
            found = true;
            firstIndex = i;
        }
        i++;
    }
    for(int j = firstIndex;j<graph->getVertices();j++){
        cout<<bestWay[j]<<" ";
    }
    for (int k = 0; k < firstIndex; ++k) {
        cout<<bestWay[k]<<" ";
    }
    cout<<"Cost: "<<bestCost<<endl;
    cout<<"bruteforce time: "<<algorithmTime<<endl;

}