//
// Created by ado on 30.10.15.
//

#ifndef TSP_ALGORITHMS_H
#define TSP_ALGORITHMS_H


class Algorithms {
public:

    double getTime();
    virtual void runAlgorithm() =0;
    virtual void printResult() =0;
    int getBestCost();
    int* getBestWay();
protected:
    double algorithmTime;
    int* bestWay;
    int bestCost;
};


#endif //TSP_ALGORITHMS_H
