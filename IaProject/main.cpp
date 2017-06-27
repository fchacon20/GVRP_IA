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
using namespace std;

//const int N_CITIES = 25;

struct City{
    double latitude;
    double longitude;
    string name;
    string type;
};

double evaluation(vector<int> sol, vector<vector<double> > d){
    double dist = 0;
    int prev = 0;

    for(vector<int>::iterator it = sol.begin() + 1; it != sol.end(); ++it) {
        dist+= d[prev][*it];
        prev = *it;
    }
    return(dist);
}

/*string initialGuess(City c[N_CITIES], double d[N_CITIES][N_CITIES], double p[4],
                        double Q, double r, double TL, double v) {
    int solution[N_CITIES];
    for (int i = 0; i < N_CITIES; ++i)
        solution[i] = i;
    bool feasibleSolution = false;
    int random, random2;
    int prev = solution[0];
    double totalDistance = 0;

    srand ((unsigned int) time(NULL));

    while (!feasibleSolution){
        //cout << endl << "new iteration" << endl;
        totalDistance = 0;
        Q = p[0];
        for (int i = 1; i < N_CITIES; ++i) {
            Q -= d[solution[prev]][solution[i]];
            //cout << "Q despues del viaje " << Q << endl;
            if (Q < 0) {
                //cout << "Inviable: Sin combustible" << endl;
                random  = rand() % 25;
                random2 = rand() % 25;
                while(random == random2)
                    random2 = rand() % 25;
                swap(solution[random],solution[random2]);
                feasibleSolution = false;
                break;
            }
            feasibleSolution = true;
            if (c[prev].type != "c")
                Q = p[0];

            totalDistance += d[solution[prev]][solution[i]];
            prev = i;
        }
    }

    cout << totalDistance << endl;
    return("");
}*/

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
    for (int i = 0; i < N_CITIES-1; ++i)
        cout << solution[i] << "-";
    cout << solution[N_CITIES-1] << endl;

    bool feasibleSolution = false;
    int prev = solution[0];
    double totalDistance = 0;
    vector<int> vectorSolution;
    vectorSolution.push_back(solution[0]);
    double goBack = 0;
    int time = 0;

    while (!feasibleSolution){
        totalDistance = 0;
        Q = p[0];

        for (int i = 1; i < N_CITIES; ++i) {
            goBack = d[solution[prev]][solution[0]];
            Q -= r*d[solution[prev]][solution[i]];

            if ((Q < goBack)&&(goBack != 9999)) {
                cout << "Inviable: Sin combustible en ciudad " << c[solution[prev]].name << endl;
                vectorSolution.push_back(solution[0]);
                totalDistance += d[solution[prev]][solution[0]];
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

            if(time >= 60*TL){
                cout << "Limite de tiempo alcanzado" << endl;
                vectorSolution.push_back(solution[0]);
                totalDistance += d[solution[prev]][solution[0]];
                prev = solution[0];
                i--;
                time = 0;
                Q = p[0];
                continue;
            }

            totalDistance += d[solution[prev]][solution[i]];
            vectorSolution.push_back(solution[i]);
            prev = i;
        }
    }

    cout << "Tiempo total: " << time << endl;

    for(vector<int>::iterator it = vectorSolution.begin(); it != vectorSolution.end(); ++it) {
        cout << *it << "-" ;
    }
    cout << endl;
    cout << "Distancia de solucion inicial: " << totalDistance << endl;
    return(vectorSolution);
}

bool isFeasible(vector<City> c, vector<int> sol, vector<vector<double> > d, double p[4]){
    double Q = p[0];
    double r = p[1];
    double TL = p[2];
    //double v = p[3];
    int N_CITIES = (int)c.size();

    int prev = sol[0];
    int time = 0;

    for (int i = 1; i < N_CITIES; ++i) {
        Q -= r*d[sol[prev]][sol[i]];

        if (Q < 0) {
//            cout << "Inviable: Sin combustible en ciudad " << c[sol[prev]].name << endl;
            return false;
        }

        if (c[prev].type != "c") {
            Q = p[0];
            time += 15;
        }else
            time += 30;

        if(time >= 60*TL){
//            cout << "Limite de tiempo alcanzado" << endl;
            return false;
        }
        prev = i;
    }

    return true;
}

double toRadians(double degrees){
    return((degrees * M_1_PI)/180);
}

double haversine(double lat1, double lat2, double lon1, double lon2) {
    double radiusOfEarth = 4182.44949; // miles, 6371km;
    double dLat = toRadians(lat2 - lat1);
    double dLon = toRadians(lon2 - lon1);
    double a = sin(dLat / 2) * sin(dLat / 2) + cos(toRadians(lat1)) *
                                                         cos(toRadians(lat2)) * sin(dLon / 2) *
                                                         sin(dLon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    double distance = radiusOfEarth * c;

    return distance;
}

int main(int argc, char** argv) {

    string temp(argv[1]);
    int iterations = atoi(argv[2]);
    //City cities[N_CITIES];
    vector<City> cities;
    string line;
    ifstream file(temp);
    vector<string> text;
    int index = 0;
    vector<int> solution;

    //double distances[N_CITIES][N_CITIES];
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

    //Matriz de distancia
    /*for (int i = 0; i < N_CITIES; ++i) {
        for (int j = 0; j < N_CITIES; ++j) {
            cout << distances[i][j] << " ";
        }
        cout << endl;
    }*/
    //cout << distances[0][0] << endl << endl;

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

    srand ((unsigned int) time(NULL));

    vector<int> bestSolution = solution;
    double bestDistance = distPrev;
    int temperature = 400;
    double dif = 0;
    double prob = 0;
    double ran = 0;
    int incremental = 0;

    while(incremental < iterations) {
        for (unsigned int i = 1; i < solution.size() - 1; ++i) {
            for (unsigned int j = 1; j < solution.size() - 1; ++j) {
                if (i != j) {
                    iter_swap(solution.begin() + i, solution.begin() + j);

                    //Si no es factible, se disolve el swap
                    if (!isFeasible(cities, solution, distances, parameters)) {
                        iter_swap(solution.begin() + i, solution.begin() + j);
                    } else {
                        newDist = evaluation(solution, distances);

                        if (bestDistance > newDist) {
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
                //cout << newDist << endl;
            }
            if (incremental > iterations) break;
        }
    }

    for(vector<int>::iterator it = bestSolution.begin(); it != bestSolution.end(); ++it) {
        cout << *it << "-" ;
    }
    cout << endl << "la mejor solucion es: " << evaluation(bestSolution, distances) << endl;

    return 0;
}




