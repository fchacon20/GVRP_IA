#include <cmath>
#include <random>
#include <utility>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <cstring>
#include <math.h>
#include <time.h>
#include "haversine.cpp"
#include "solutionUtilities.cpp"
#include "city.h"
using namespace std;

vector<int> initialGuess(vector<City> c, vector<vector<double> > d, double p[4],
                        double Q, double r, double TL, double v) {
    int N_CITIES = (int)c.size();
    int solution[N_CITIES];
    double minimun = 0;
    vector<int> ready;
    int next = 0;

    bool feasible = false;
    ready.push_back(0);
    solution[0] = 0;
    for (int i = 1; i < N_CITIES-1; ++i) {
        vector<double> auxV;
        for (int j = 0; j < N_CITIES; ++j)
            auxV.push_back(d[solution[i-1]][j]);

        sort(auxV.begin(), auxV.end());

        while(!feasible) {
            minimun = auxV[next];
            if (minimun == 9999) {
                solution[i] = solution[i-1];
                break;
            }
            for (int j = 0; j < N_CITIES; ++j) {
                if (d[solution[i - 1]][j] == minimun) {
                    if(find(ready.begin(),ready.end(),j) != ready.end()) {
                        next++;
                        feasible = false;
                        break;
                    }
                    ready.push_back(j);
                    solution[i] = j;
                }
                feasible = true;
            }
        }
        feasible = false;
        next = 0;
    }

    solution[N_CITIES-1] = 0;

    bool feasibleSolution = false;
    int prev = solution[0];
    vector<int> vectorSolution;
    vectorSolution.push_back(solution[0]);
    double goBack = 0;
    int time = 0;

    while (!feasibleSolution){
        Q = p[0];

        for (int i = 1; i < N_CITIES; ++i) {
            if(solution[prev] == solution[i])
                continue;

            goBack = d[solution[prev]][solution[0]];
            Q -= r*d[solution[prev]][solution[i]];

            //Sin combustible
            if ((Q < goBack)&&(goBack != 9999)) {
                vectorSolution.push_back(solution[0]);
                prev = solution[0];
                i--;
                time = 0;
                Q = p[0];
                continue;
            }

            feasibleSolution = true;
            if (c[prev].type != "c") {
                Q = p[0];
                time += 15;
            }else
                time += 30;

            time += d[solution[prev]][solution[i]]/v;

            //Tiempo excedido
            if(time >= 60*TL){
                vectorSolution.push_back(solution[0]);
                prev = solution[0];
                i--;
                time = 0;
                Q = p[0];
                continue;
            }

            vectorSolution.push_back(solution[i]);
            prev = i;
        }
        vectorSolution.push_back(0);

    }

    printVector(vectorSolution, c);
    cout << "Distancia de solucion inicial: " << evaluation(vectorSolution, d) << endl << endl;
    return(vectorSolution);
}

int main(int argc, char** argv) {

    clock_t begin = clock();
    string temp(argv[1]);
    int iterations = atoi(argv[2]);
    vector<City> cities;
    string line;
    ifstream file(temp);
    vector<string> text;
    int index = 0;
    vector<int> solution;
    vector<vector<double> > distances;

    getline(file, line);
    while ((getline(file, line)) && (line.length() > 1)) {
        vector<string> result;
        istringstream iss(line);
        for (string s; iss >> s;)
            result.push_back(s);

        cities.push_back(City());
        cities[index].name = result[0];
        cities[index].type = result[1];
        cities[index].longitude = stod(result[2]);
        cities[index].latitude = stod(result[3]);
        index++;
    }
    int N_CITIES = index;
    double aux;


    for (int i = 0; i < N_CITIES; ++i) {
        vector<double > row; //Empty row
        for (int j = 0; j < N_CITIES; ++j) {
            row.push_back(0);
        }
        distances.push_back(row);
    }

    for (int i = 0; i < N_CITIES; ++i) {
        distances[i][i] = 9999;
        for (int j = 0; j < i; ++j) {
            aux = haversine(cities[i].latitude,cities[j].latitude, cities[i].longitude, cities[j].longitude);
            if (aux == 0)
                aux = 9999;
            distances[i][j] = aux;
            distances[j][i] = aux;
        }
    }

    //Parametros
    double parameters[4];
    index = 0;
    while(getline(file,line)) {
        if (line.length() <= 1)
            continue;
        size_t found1 = line.find("/");
        size_t found2 = line.find("/",found1+1,1);
        if (found2!=std::string::npos)
            line = line.substr(found1+1,found2-found1-1);

        parameters[index] = stod(line);
        index++;
    }

    double Q   = parameters[0];//gallons
    double r   = parameters[1];// gallons per mile
    double TL  = parameters[2];// hours
    double v   = parameters[3];//miles per hours
    v = v/60; //miles per minutes

    solution = initialGuess(cities, distances, parameters, Q, r, TL, v);
    double distPrev = evaluation(solution, distances);
    double newDist = 0;
    vector<int> firstSolution = solution;

    srand ((unsigned int) time(NULL));

    vector<int> bestSolution = solution;
    double bestDistance = distPrev;
    int bestIteration = 0;
    int temperature = 400;
    double dif = 0;
    double prob = 0;
    double ran = 0;
    int incremental = 0;
    int otherRandom = 0;

    while(incremental < iterations) {
        for (unsigned int i = 1; i < solution.size() - 1; ++i) {
            for (unsigned int j = 1; j < solution.size() - 1; ++j) {
                if (i != j) {
                    iter_swap(solution.begin() + i, solution.begin() + j);

                    //Si no es factible, se disuelve el swap
                    if (!isFeasible(cities, solution, distances, parameters)) {
                        iter_swap(solution.begin() + i, solution.begin() + j);
                    } else {
                        newDist = evaluation(solution, distances);

                        if (bestDistance > newDist) {
                            bestIteration = incremental;
                            bestSolution = solution;
                            bestDistance = evaluation(bestSolution, distances);
                        }

                        //Si la nueva solucion es mejor que la anterior, se deja pasar
                        if (distPrev >= newDist) {
                            distPrev = newDist;
                        } else if (temperature > 0) {
                            dif = newDist - distPrev;
                            prob = exp(dif / temperature);
                            ran = ((double) rand() / (RAND_MAX));

                            //Si la exponencial es menor que el numero aleatorio, no
                            //se deja pasar, es decir, se disuelve el swap
                            if (ran > prob) {
                                iter_swap(solution.begin() + i, solution.begin() + j);
                            } else {
                                distPrev = newDist;
                                temperature -= 1;
                            }
                        }
                    }
                }
                incremental++;
                if (incremental > iterations) break;
            }
            if (incremental > iterations) break;
        }

        if (solution == firstSolution) {
            otherRandom = (int) ((rand() % (solution.size() - 2)) + 2);
            solution.push_back(0);
            iter_swap(solution.begin() + otherRandom , solution.end() - 2);
            firstSolution = solution;
        }
    }

    cout << endl;
    printVector(bestSolution, cities);
    cout << endl << "La mejor solucion es: " << evaluation(bestSolution, distances) << endl;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    cout << "Solucion encontrada en la iteracion: " << bestIteration << endl;
    cout << "Tiempo de ejecucion: " << time_spent << endl;

    return 0;
}




