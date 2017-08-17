#ifndef IAPROJECT_SOLUTIONUTILITIES_H
#define IAPROJECT_SOLUTIONUTILITIES_H
#include "city.h"

double evaluation(std::vector<int> sol, std::vector<std::vector<double> > d);
void printVector(std::vector<int> sol);
bool isFeasible(std::vector<City> c, std::vector<int> sol, std::vector<std::vector<double> > d, double p[4]);

#endif //IAPROJECT_SOLUTIONUTILITIES_H
