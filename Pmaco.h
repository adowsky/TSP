#ifndef TSP_PMACO_H
#define TSP_PMACO_H
#include "Graph.h"
#include "Algorithms.h"
class Pmaco : public Algorithms {
public:
    Pmaco(Graph *g);

    void runAlgorithm();

    void runAlgorithm(int travelsAmount);

    void printResult();

private:
    Graph *graph;
    int **normalAntsPopulation;
    int **hybridAntsPopulation;
    int mutatedPopulationAmount;
    double **pheromones;
    double **delta;
    bool *visited;
    int localUpdateLimit;

    const int ALPHA_FACTOR = 3;
    const int BETA_FACTOR = 2;
    const double EVAPORATION_FACTOR = 0.01;
    const double HYBRID_POPULATION_FACTOR_MIN = 0.05;
    const double HYBRID_POPULATION_FACTOR_MAX = 0.15;
    const double PHEROMONE_INITIAL_VALUE = 5.0;
    const double INCREASE_FACTOR = 1.2;

    void initialize();
    bool *createVisitedTab();
    void generateNormalAntsRoutes(int iteration);
    int getVertexWithBestProbability(int vertex);
    void singlePheromoneUpdate(int iteration, int vertexFrom, int vertexTo);
    void localPheromonesUpdate(int vertexFrom, int vertexTo);
    int nearestNeighbourLength(int vertex);
    void deltaUpdate(int from, int to);
    void clearVisitedTab();
    void setBestWay(int *way, int cost);
    void deltaUpdate(int* way);
    void globalPheromonesUpdate();
    void mutation(int currentPopulationAmount);
};

#endif //TSP_PMACO_H
