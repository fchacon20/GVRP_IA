#include <math.h>
#include "haversine.h"

//Funcion que transforma grados en radianes
double toRadians(double degrees){
    return((degrees * M_1_PI)/180);
}

//Funcion que retorna la distancia entre dos nodos dadas sus latitudes y longitudes
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
