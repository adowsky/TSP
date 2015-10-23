#include <iostream>
#include<algorithm>
#include<fstream>
#include <cmath>


using namespace std;
/**
 * Graph representation
 */
class Graph{
public:
    Graph(string s);
    Graph(bool aut);
    void printGraph();
    void bruteForceAlgorithm();
    double getBruteForceTime();
    void antColonyAlgorithm(int antsAmount,int travelAmount);
    double getAntColonyAlgorithmTime(int antsAmount, int travelsAmount);
    void geneticAlgorithm();
    double getGeneticAlgorithmTime();
    void dynamicProgrammingAlgorithm(); // not a heuristic
    double getDynamicProgrammingAlgorithm();

private:
    int vertices;
    int** graph;
    int firstVertex;
};
Graph::Graph(string s) {
    fstream file;
    file.open(s,ios::in);
    if(file.good()){
//        reading graph details
        file>>this->firstVertex;
        file>>this->vertices;
        graph=new int*[vertices];
        for (int i = 0; i < vertices; ++i) {
            graph[i]=new int[vertices];
        }
//      reading graph representation and writing it to matrix
        for (int i = 0; i < vertices; ++i) {
            for (int j = 0; j < vertices; ++j) {
                file>>graph[i][j];
            }
        }
    }
    file.close();

}
Graph::Graph( bool aut){
    cout<<"Vertices: ";
    cin>>vertices;
    graph= new int*[vertices];
    for(int i=0;i<vertices;++i)
        graph[i]=new int[vertices];
    srand(time(NULL));
    for (int j = 0; j < vertices; ++j) {
        for (int i = 0; i < vertices; ++i) {
            if(j==i)
                graph[i][j]=0;
            else if(j>i){
                graph[j][i]=graph[i][j];
            }
            else{
                graph[j][i]= rand()%(vertices*10) +1;
            }

        }
    }
    firstVertex=rand()%vertices;
}
void Graph::bruteForceAlgorithm() {
    int sequence[vertices-1];
    int i=0;
    int permutationAmount=0;
    int bestWay[vertices];
    bestWay[0]=firstVertex;
    int bestCost=0;
    for (int j = 0,i=0; j < vertices-1; ++j,++i) {
     if(i==firstVertex){
         j--;
         continue;
     }
        sequence[j]=i;
    }
//    loop which checks every next permutation of possible result
    do{
        permutationAmount++;
        int cost=graph[firstVertex][sequence[0]];
        if(cost==0) continue;
        for (int i = 0; i < vertices-2; ++i) {
            if(graph[sequence[i]][sequence[i+1]]==0){
                cost=0;
                break;
            }
            cost+=graph[sequence[i]][sequence[i+1]];
        }
        cost+=graph[sequence[vertices-2]][firstVertex];
        if(cost<=bestCost||bestCost==0) {
            bestCost = cost;
//            rewriting sequence
            for (int j = 0; j < vertices-1; ++j) {
                bestWay[j+1]=sequence[j];
            }
        }
    }while(next_permutation(sequence,sequence+vertices-1));
    cout<<"Checked "<<permutationAmount<<" permutations"<<endl;
    for (int k = 0; k < vertices; ++k) {
        cout<<bestWay[k];
    }
    cout<<endl<<"Best cost: "<<bestCost;
}
double Graph::getBruteForceTime() {
    time_t time=clock();
    bruteForceAlgorithm();
    return (clock()-time)/(double)CLOCKS_PER_SEC;
}

/**
 * Algorithm random path depending on factors between every two vertices that mark best solution.
 * @parm antsAmount number of algorithm's iterations
 */
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
double Graph::getAntColonyAlgorithmTime(int antsAmount, int travelsAmount) {
    clock_t time = clock();
    antColonyAlgorithm(antsAmount,travelsAmount);
    return (clock()-time)/(double)(CLOCKS_PER_SEC);
    
}
int main() {
    //Graph g = Graph("graph.txt");
    Graph g=Graph(false);
   // g.printGraph();
    cout<<endl<<"Brute force time: "<<g.getBruteForceTime()<<endl;
    cout<<"Ants Colony Time: "<<g.getAntColonyAlgorithmTime(3,1000);
    return 0;
}

void Graph::printGraph() {
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            cout<<graph[i][j]<<" ";
        }
        cout<<endl;
    }

}
