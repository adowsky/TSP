#include "bruteForce.h"
#include "Pmaco.h"


using namespace std;

int main() {
    Graph* g = new Graph();
    g->printGraph();
    bruteForce* alg = new bruteForce(g);
    //alg->runAlgorithm();
    //alg->printResult();
    //int best = alg->getBestCost();
    Algorithms *alg2 = new Pmaco(g);
    ((Pmaco*)alg2)->classicAnts();
    alg2->printResult();
    alg2->runAlgorithm();
    alg2->printResult();
    return 0;
}
