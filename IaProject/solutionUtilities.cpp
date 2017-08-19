#include <vector>
#include <iostream>
#include "solutionUtilities.h"
#include "city.h"
using namespace std;

//Evaluacion en la funcion objetivo de la solucion propuesta
double evaluation(vector<int> sol, vector<vector<double> > d){
    double dist = 0;
    int prev = 0;
    for(vector<int>::iterator it = sol.begin() + 1; it != sol.end() - 1; ++it) {
        dist+= d[prev][*it];
        prev = *it;
    }
    return(dist);
}

//Imprimir la solucion de manera legible
void printVector(vector<int> sol, vector<City> c){
    for(vector<int>::iterator it = sol.begin(); it != sol.end() - 2; ++it)
        cout << c[*it].name << "-";
    cout << c[*(sol.end()-1)].name << endl;
}

//Comprueba si la solucion es factible, es decir, si se satisfacen las restricciones
//de combustible y tiempo
bool isFeasible(vector<City> c, vector<int> sol, vector<vector<double> > d, double p[4]){
    double Q = p[0];
    double r = p[1];
    double TL = p[2];
    double v = p[3]; //miles per hour
    v = v/60; //miles per minute
    int N_CITIES = (int)c.size();

    int prev = sol[0];
    int time = 0;

    for (int i = 1; i < N_CITIES; ++i) {
        Q -= r * d[sol[prev]][sol[i]];

        if (Q < 0) {
            return false;
        }

        if (c[sol[prev]].type != "c") {
            Q = p[0];
            time += 15;
        } else {
            time += 30;
        }

        if(c[sol[prev]].type == "d"){
            Q = p[0];
            time = 0;
        }

        time += d[sol[prev]][sol[i]]/v;

        if(time >= 60*TL) {
            return false;
        }

        prev = i;
    }

    return true;
}