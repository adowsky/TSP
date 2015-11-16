#include "bruteForce.h"
#include "Pmaco.h"


using namespace std;
/*
void Graph::antColonyAlgorithm(int antsAmount, int travelAmount) {
    int alphaValue = 3; // how important pheromone is
    int betaValue =2;
    double minimumPheromoneValue= 0.01;
    double initialValueOfPheromone =20.0;
    double evaporationFactor =0.01; // x: (0,1>
    double pheromoneIncreaseFactor =2;
    double maxPheromoneValue=10000.0;
    int rouletteFactor =90;
    srand(time(NULL));
    double **pheromones = new double*[vertices];
    for (int m = 0; m < vertices; ++m) {
        pheromones[m]=new double[vertices];
    }
    double** deltaPheromones =new double*[vertices];
    for (int l1 = 0; l1 < vertices; ++l1) {
        deltaPheromones[l1]=new double[vertices];
    }
    for (int n = 0; n < vertices; ++n) {
        for (int i = 0; i < vertices; ++i) {
            if(n!=i) {
                pheromones[n][i] = initialValueOfPheromone;
                deltaPheromones[n][i]=0;
            }
        }
    }
    int bestCost=0;
    int* bestWay = new int[vertices];
    for (int i = 0; i < travelAmount; ++i) { //travels
        int** ants = new int*[antsAmount];
        for (int j = 0; j < antsAmount; ++j) {
            ants[j]=new int[vertices];
        }
        for (int j = 0; j < antsAmount; ++j) { //every ant
            //ants[j][0]=firstVertex;
            ants[j][0]= rand()%vertices;
            bool* visited =new bool[vertices];
            for (int i1 = 0; i1 < vertices; ++i1) {
                visited[i1]=false;
            }
            visited[ants[j][0]]=true;
            int cost=0;
                for (int k = 1; k < vertices; ++k) { //every vertex
                    int nextVertex = -1;
                    double sum = 0;
                    for (int l = 0; l < vertices; ++l) {//counting sum of unvisited vertices
                        if (!visited[l] && l != k)
                            sum +=(pow(pheromones[ants[j][k-1]][l], alphaValue) * pow(((double)1 / graph[ants[j][k-1]][l]), betaValue));
                    }//sum
                    int bestProb = 0;
                    for (int m = 0; m < vertices; ++m) { // finding best way
                        if (visited[m]) continue;
                        int roulette = rand() % 100 + 1;
                        double probability=0.0;
                        if (roulette < rouletteFactor) {
                            probability = pow(pheromones[ants[j][k - 1]][m], alphaValue) *
                                                 pow(((double) 1 / graph[ants[j][k - 1]][m]), betaValue) / sum;

                        }
                        if (probability >= bestProb) {
                            bestProb = probability;
                            nextVertex = m;
                        }
                    }//best way
                    //updating ant
                    cost+=graph[ants[j][k-1]][nextVertex];
                    ants[j][k]=nextVertex;
                    visited[nextVertex]=true;
                } //every vertex
            cost+=graph[ants[j][vertices-1]][ants[j][0]];
            if(cost<=bestCost||bestCost==0){ // verifying cost
                bestCost=cost;
                for (int k = 0; k < vertices; ++k) {
                  //  cout<<ants[j][k];
                    bestWay[k]=ants[j][k];
                }
                //cout<<endl;
            }

        }//every ant
        for (int j = 0; j < antsAmount; ++j) { // calculating sum of length of every's ant route
            double length=0.0;
            for (int l = 0; l < vertices-1; ++l) {
                length+=graph[ants[j][l]][ants[j][l+1]];
            }
            length+=graph[ants[j][vertices-1]][ants[j][0]];
            for (int k = 0; k < vertices-1; ++k) {
                deltaPheromones[ants[j][k]][ants[j][k+1]]+=pheromoneIncreaseFactor/length;
            }
            deltaPheromones[ants[j][vertices-1]][ants[j][0]]+=pheromoneIncreaseFactor/length;
        }
        for (int n = 0; n < vertices; ++n) { //pheromones update
            for (int k = 0; k < vertices; ++k) {
                if (n != k) {
                double val = (1 - evaporationFactor) * pheromones[n][k] + deltaPheromones[n][k];
                    deltaPheromones[n][k]=0;
                if (val < minimumPheromoneValue)
                    pheromones[n][k] = minimumPheromoneValue;
                else if (val > maxPheromoneValue)
                    pheromones[n][k] = maxPheromoneValue;
                else
                    pheromones[n][k] = val;
                }
            }
        }

    }//travels
    //printing best of
    cout<<"Best way: ";
    for (int k1 = 0; k1 < vertices; ++k1) {
        cout<<bestWay[k1];
    }
    cout<<endl;
    cout<<"cost: "<<bestCost<<endl;
}//end of AntsColonyAlgorithm
*/

int main() {
    Graph* g = new Graph();
    g->printGraph();
    bruteForce* alg = new bruteForce(g);
    alg->runAlgorithm();
    alg->printResult();
    int best = alg->getBestCost();
    Algorithms *alg2 = new Pmaco(g);
    alg2->runAlgorithm();
    alg2->printResult();
    cout<<"Przybliżenie wynosi: "<<alg2->getBestCost()/(double)best;
    return 0;
}
