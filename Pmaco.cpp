#include "Pmaco.h"
#include <cmath>

Pmaco::Pmaco(Graph *g) {
    graph=g;
}

void Pmaco::initialize() {
    int vertices = graph->getVertices();
    mutatedPopulationAmount = vertices * HYBRID_POPULATION_FACTOR_MAX;
    if( (mutatedPopulationAmount % 2) == 1 )
        mutatedPopulationAmount++;
    createVisitedTab();
    bestCost = 0;
    localUpdateLimit = vertices / 5;
    normalAntsPopulation = new int *[vertices];
    hybridAntsPopulation = new int *[mutatedPopulationAmount];
    pheromones = new double *[vertices];
    delta = new double *[vertices];
    bestWay = new int[vertices];
    for (int i = 0; i < vertices; ++i) {
        normalAntsPopulation[i] = new int[vertices];
        pheromones[i] = new double[vertices];
        delta[i] = new double[vertices];
        bestWay[i] = -1;
        for (int j = 0; j < vertices; ++j) {
            normalAntsPopulation[i][j] = -1;
            pheromones[i][j] = PHEROMONE_INITIAL_VALUE;
            delta[i][j] = 0;
        }
    }
    for (int k = 0; k < mutatedPopulationAmount; ++k) {
        hybridAntsPopulation[k] = new int[vertices];
        for (int i = 0; i < vertices; ++i) {
            hybridAntsPopulation[k][i] = -1;
        }
    }

}

int Pmaco::getVertexWithBestProbability(int vertex) {
    double bestProb = 0;
    int possibleVertex = -1;
    int vertices = graph->getVertices();
    double probabilityContainer[vertices];
    double probSum = 0;
    int **routesGraph = graph->getGrpah();
    for (int i = 0; i < vertices; ++i) {
        if (visited[i]) {
            probabilityContainer[i] = 0;
            continue;
        }
        // Computing value using specific equation
        probabilityContainer[i] = pow(pheromones[vertex][i], ALPHA_FACTOR) *
                                  pow(((double) 1 / routesGraph[vertex][i]), BETA_FACTOR);
        probSum += probabilityContainer[i];
        if (probabilityContainer[i] > bestProb) {
            bestProb = probabilityContainer[i];
            possibleVertex = i;
        }
    }
    if (possibleVertex == -1)
        cout << "Error while computing probability!" << endl;
    return possibleVertex;
}

bool *Pmaco::createVisitedTab() {
    visited = new bool[graph->getVertices()];
    for (int i = 0; i < graph->getVertices(); ++i) {
        visited[i] = false;
    }
    return visited;
}

void Pmaco::clearVisitedTab() {
    for (int i = 0; i < graph->getVertices(); ++i) {
        visited[i] = false;
    }
}

void Pmaco::singlePheromoneUpdate(int iteration, int vertexFrom, int vertexTo) {
    if (iteration < localUpdateLimit)
        localPheromonesUpdate(vertexFrom, vertexTo);
    else
        deltaUpdate(vertexFrom, vertexTo);

}

void Pmaco::localPheromonesUpdate(int vertexFrom, int vertexTo) {
    pheromones[vertexFrom][vertexTo] = (1 - EVAPORATION_FACTOR) * pheromones[vertexFrom][vertexFrom] +
                                       EVAPORATION_FACTOR *
                                       (1.0 / (graph->getVertices() * nearestNeighbourLength(vertexFrom)));
}

int Pmaco::nearestNeighbourLength(int vertex) {
    int cost = 0;
    int **routeGraph = graph->getGrpah();
    int vertices = graph->getVertices();
    bool visit[vertices];
    for (int k = 0; k < vertices; ++k) {
        visit[k] = false;
    }

    for (int i = 0; i < vertices; ++i) {
        visit[vertex] = true;
        int bestNext = -1;
        int lowestCost = 0;
        for (int j = 0; j < vertices; ++j) {
            if (visit[j]) continue;
            if (routeGraph[vertex][j] < lowestCost || lowestCost == 0) {
                bestNext = j;
                lowestCost = routeGraph[vertex][j];
            }
        }
        cost += lowestCost;
        visit[bestNext] = true;
        vertex = bestNext;
    }
    return cost;
}

void Pmaco::deltaUpdate(int from, int to) {
    delta[from][to] += INCREASE_FACTOR;
}

void Pmaco::setBestWay(int *way, int cost) {
    for (int i = 0; i < graph->getVertices(); ++i) {
        bestWay[i] = way[i];
    }
    bestCost = cost;
}

void Pmaco::generateNormalAntsRoutes(int iteration) {
    int **routesGraph = graph->getGrpah();
    for (int i = 0; i < graph->getVertices(); ++i) {
        int cost = 0;
        normalAntsPopulation[i][0] = i;
        visited[i] = true;
        for (int j = 1; j < graph->getVertices(); ++j) { //every vertex
            normalAntsPopulation[i][j] = getVertexWithBestProbability(j - 1);
            visited[normalAntsPopulation[i][j]] = true;
            cost += routesGraph[normalAntsPopulation[i][j - 1]][normalAntsPopulation[i][j]];
            singlePheromoneUpdate(iteration, normalAntsPopulation[i][j - 1], normalAntsPopulation[i][j]);
        }
        cost += routesGraph[normalAntsPopulation[i][graph->getVertices()-1]][normalAntsPopulation[i][0]];
        if (cost < bestCost || bestCost == 0) {
            setBestWay(normalAntsPopulation[i], cost);
        }
        clearVisitedTab();
    }
}

void Pmaco::globalPheromonesUpdate() {
    int vertices = graph->getVertices();
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            double val = pheromones[i][j] * (1 - EVAPORATION_FACTOR) + delta[i][j];
            delta[i][j] = 0;
            if (val < EVAPORATION_FACTOR)
                pheromones[i][j] = EVAPORATION_FACTOR;
            else
                pheromones[i][j] = val;

        }
    }
}

void Pmaco::runAlgorithm(int travelsAmount) {
    clock_t timer = clock();
    initialize();
    if((HYBRID_POPULATION_FACTOR_MIN + HYBRID_POPULATION_FACTOR_MAX-HYBRID_POPULATION_FACTOR_MIN)*graph->getVertices() > 1)
        cout<<"Hybrid population is used in algorithm!"<<endl;
    for (int i = 0; i < travelsAmount; ++i) {
        generateNormalAntsRoutes(i);
        double currentHybridFactor = HYBRID_POPULATION_FACTOR_MIN + (i/(double)travelsAmount)*
                                                                            (HYBRID_POPULATION_FACTOR_MAX-HYBRID_POPULATION_FACTOR_MIN);
        if((int)(graph->getVertices()*currentHybridFactor) % 2 ==1)
            currentHybridFactor--;
        if(graph->getVertices() * currentHybridFactor >2)
            mutation(currentHybridFactor*graph->getVertices());
            //generateMutatedAntsRoutes();
        globalPheromonesUpdate();
    }
    algorithmTime = (clock() - timer) / (double) CLOCKS_PER_SEC;
}
void Pmaco::classicAnts() {
    clock_t timer = clock();
    initialize();
    for (int i = 0; i < graph->getVertices()*graph->getVertices(); ++i) {
        generateNormalAntsRoutes(i);
        globalPheromonesUpdate();
    }
    algorithmTime = (clock() - timer) / (double) CLOCKS_PER_SEC;

}

void Pmaco::runAlgorithm() {
    runAlgorithm(graph->getVertices()*graph->getVertices());

}

void Pmaco::deltaUpdate(int *route) {
    int vertices = graph->getVertices();
    for (int i = 0; i <vertices; ++i) {
        delta[route[i]][route[(i + 1) % vertices]] += INCREASE_FACTOR;
    }
    }
    void Pmaco::mutation(int currentPopulationAmount){

        for (int i = 0; i < currentPopulationAmount; i+=2) { // for every ant
            if(currentPopulationAmount>=4){
                int z=0;
            }
            int firstCost = 0;
            int secondCost = 0;
            int vertices = graph->getVertices();
            int beginVertex = rand()%vertices;
            int endVertex = rand()%vertices;
            while(endVertex == beginVertex)
                endVertex = rand()%vertices;
            if(endVertex<beginVertex){
                int tmp = endVertex;
                endVertex = beginVertex;
                beginVertex = tmp;
            }
            int firstParent = rand()%vertices;
            int secondParent = rand()%vertices;
            while(firstParent == secondParent)
                secondParent = rand()%vertices;
            for (int j = beginVertex; j <= endVertex ; ++j) {    //replacing specified elements of route
                hybridAntsPopulation[i][j] = normalAntsPopulation[secondParent][j];
                hybridAntsPopulation[i+1][j] = normalAntsPopulation[secondParent][j];
            }
            //finding fitting elements
            for (int k = 0; k < beginVertex; ++k) {
                bool firstFitting = true;
                bool secondFitting = true;
                for (int j = beginVertex; j <= endVertex; ++j) {
                    if(normalAntsPopulation[firstParent][k] == hybridAntsPopulation[i][j]){
                        firstFitting = false;
                        if(!secondFitting)
                            break;
                    }
                    if(normalAntsPopulation[secondParent][k] == hybridAntsPopulation[i+1][j]){
                        secondFitting = false;
                        if (!firstFitting)
                            break;
                    }
                }
                if(firstFitting){
                    hybridAntsPopulation[i][k] = normalAntsPopulation[firstParent][k];
                }
                else
                    hybridAntsPopulation[i][k] = -1;
                if(secondFitting)
                    hybridAntsPopulation[i+1][k] = normalAntsPopulation[secondParent][k];
                else
                    hybridAntsPopulation[i+1][k] = -1;
            }
            for (int k = endVertex+1; k < vertices; ++k) {
                bool firstFitting = true;
                bool secondFitting = true;
                for (int j = beginVertex; j <= endVertex; ++j) {
                    if(normalAntsPopulation[firstParent][k] == hybridAntsPopulation[i][j]){
                        firstFitting = false;
                        if(!secondFitting)
                            break;
                    }
                    if(normalAntsPopulation[secondParent][k] == hybridAntsPopulation[i+1][j]){
                        secondFitting = false;
                        if (!firstFitting)
                            break;
                    }
                }
                if(firstFitting){
                    hybridAntsPopulation[i][k] = normalAntsPopulation[firstParent][k];
                }
                else
                    hybridAntsPopulation[i][k] = -1;
                if(secondFitting)
                    hybridAntsPopulation[i+1][k] = normalAntsPopulation[secondParent][k];
                else
                    hybridAntsPopulation[i+1][k] = -1;
            }
            //filling missing elements
            for (int l = 0; l < beginVertex; ++l) {
                if(hybridAntsPopulation[i][l] == -1){
                    for (int j = 0; j < vertices; ++j) {
                        bool isIn = false;
                        for (int k = 0; k < vertices; ++k) {
                            if(hybridAntsPopulation[i][k] == normalAntsPopulation[secondParent][j]){
                                isIn = true;
                                break;
                            }
                        }
                        if(!isIn){
                            hybridAntsPopulation[i][l] = normalAntsPopulation[secondParent][j];
                            break;
                        }
                    }
                    for (int j = 0; j < vertices; ++j) {
                        bool isIn = false;
                        for (int k = 0; k < vertices; ++k) {
                            if(hybridAntsPopulation[i+1][k] == normalAntsPopulation[firstParent][j]){
                                isIn = true;
                                break;
                            }
                        }
                        if(!isIn){
                            hybridAntsPopulation[i+1][l] = normalAntsPopulation[firstParent][j];
                            break;
                        }
                    }
                }
            }
            for (int l = endVertex+1; l < vertices; ++l) {
                if(hybridAntsPopulation[i][l] == -1){
                    for (int j = 0; j < vertices; ++j) {
                        bool isIn = false;
                        for (int k = 0; k < vertices; ++k) {
                            if(hybridAntsPopulation[i][k] == normalAntsPopulation[secondParent][j]){
                                isIn = true;
                                break;
                            }
                        }
                        if(!isIn){
                            hybridAntsPopulation[i][l] = normalAntsPopulation[secondParent][j];
                            break;
                        }
                    }
                    for (int j = 0; j < vertices; ++j) {
                        bool isIn = false;
                        for (int k = 0; k < vertices; ++k) {
                            if(hybridAntsPopulation[i+1][k] == normalAntsPopulation[firstParent][j]){
                                isIn = true;
                                break;
                            }
                        }
                        if(!isIn){
                            hybridAntsPopulation[i+1][l] = normalAntsPopulation[firstParent][j];
                            break;
                        }
                    }
                }
            }
            int** gr = graph->getGrpah();
            for (int m = 1; m <= vertices; ++m) {
                firstCost += gr[hybridAntsPopulation[i][m-1]][hybridAntsPopulation[i][m%(vertices)]];
                secondCost += gr[hybridAntsPopulation[i+1][m-1]][hybridAntsPopulation[i+1][m%vertices]];
            }
            if(firstCost<bestCost){
                setBestWay(hybridAntsPopulation[i],firstCost);
            }
            if(secondCost<bestCost){
                setBestWay(hybridAntsPopulation[i+1],secondCost);
            }
        }
        for (int n = 0; n < currentPopulationAmount; ++n) {
            deltaUpdate(hybridAntsPopulation[n]);
        }
    }

void Pmaco::printResult() {
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
    cout<<"Pheromone-Mutation Ant Colony Optimization Time: "<<algorithmTime<<endl;

}
