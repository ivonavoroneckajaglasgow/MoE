#include "Family.h"
#include <iostream>
#include <vector>

using namespace std;
//using namespace arma;

//not what we want it to do - we need arma, which doesn't work at the moment
vector<double> Family::linkfun_vec(vector<double> x){
    vector<double> result;
    for(int i=0; i<x.size();i++)
    result[i]=this->linkfun(x[i]);
return result;
}

double Family::linkfun(double mu){

}