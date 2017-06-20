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
using namespace std;

int N_CITIES = 25;

double toRadians(double degrees){
    return((degrees * M_1_PI)/180);
}

double haversine(City city1, City city2) {
    double lat1 = city1.latitude;
    double lat2 = city2.latitude;
    double lon1 = city1.longitude;
    double lon2 = city2.longitude;
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

struct City{
    double latitude;
    double longitude;
    string name;
    string type;
};

int main(){

    City cities[N_CITIES];
    string line;
    ifstream file("test.txt");
    vector<string> text;
    int index = 0;
    vector<vector<double>> distances;

    getline(file,line);
    while((getline(file,line)) && (line != "\0")){
        vector<string> result;
        istringstream iss(line);
        for(string s; iss >> s; )
            result.push_back(s);

        cities[index].name      = result[0];
        cities[index].type      = result[1];
        cities[index].longitude = stod(result[2]);
        cities[index].latitude  = stod(result[3]);
        index++;
    }

    for (int i = 0; i < N_CITIES; ++i) {
        for (int j = 0; j < N_CITIES; ++j) {
            distances[i][j] = haversine(cities[i],cities[j]);
        }
    }

    cout << distances[0][1] << endl;
    //cout << "continuacion" << endl;
    //getline(file,line);
    //cout << line << endl;
    return 0;
}
