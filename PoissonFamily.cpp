#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "PoissonFamily.h"

using namespace std;
using namespace arma;

vec PoissonFamily::linkfun(vec mu){
    return log(mu);
}

vec PoissonFamily::linkinv(vec eta){
    return exp(eta);
}

vec PoissonFamily::var(vec mu){
    return mu;
}

double PoissonFamily::dmudeta(vec eta){
    return exp(eta).max();
}

vec PoissonFamily::density(vec y, vec eta, double lambda){
    
}