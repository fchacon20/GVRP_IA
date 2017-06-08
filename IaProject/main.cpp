//#pragma once

#include <cmath>
#include <random>
#include <utility>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>

using namespace std;

/*
// To find a status with lower energy according to the given condition
template<typename status, typename count, typename energy_function, typename temperature_function, typename next_function, typename generator>
status simulated_annealing(status i_old, count c, const energy_function& ef, const temperature_function& tf, const next_function& nf, generator& g){

    auto   e_old  = ef(i_old);

    status i_best = i_old;
    auto   e_best = e_old;

    std::uniform_real_distribution<decltype(e_old)> rf(0, 1);

    for(; c > 0; --c){
        status i_new = nf(i_old, g);
        auto   e_new = ef(i_new);

        if(e_new < e_best){
            i_best =           i_new ;
            e_best =           e_new ;
        }

        if( e_new < e_old || std::exp( (e_old - e_new) / tf(c) ) > rf(g) ){
            i_old  = std::move(i_new);
            e_old  = std::move(e_new);
        }
    }
    return(i_best);
}
class tour{
}
void SA(int T, ){
}
*/

int main(){
    string line;
    ifstream infile;
    infile.open("test.txt");

    if(infile.is_open()){
        while(getline(infile,line)){
            cout << line << endl;
            istringstream iss(line);
            vector<string> tokens{istream_iterator<string>{iss},
                                  istream_iterator<string>{}};
        }
        infile.close();
    }

    return 0;
}
