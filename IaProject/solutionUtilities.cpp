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
    int N_CITIES = (int)sol.size();

    int prev = sol[0];
    int time = 0;

    for (int i = 1; i < N_CITIES; ++i) {
        if(d[sol[prev]][sol[i]] == 9999)
            continue;

        Q -= r * d[sol[prev]][sol[i]];
        time += d[sol[prev]][sol[i]]/v;

        //cout << c[sol[i]].name  << ": " << Q << endl;
        if (Q < 0) {
            //cout << "no comb " << Q << " en " << c[sol[i]].name << endl;
            return false;
        }

        if (c[sol[i]].type == "f") {
            Q = p[0];
            time += 15;
        } else if (c[sol[i]].type == "c"){
            time += 30;
        }else{
            Q = p[0];
            time = 0;
        }

        if(time >= 60*TL) {
            //cout << "falta time: " << time << " en " << c[sol[prev]].name << endl;
            return false;
        }

        prev = i;
    }

    return true;
}