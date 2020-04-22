#define _USE_MATH_DEFINES
#define EPS 1e-9

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "GammaFamily.h"

using namespace std;
using namespace arma;

GammaFamily::GammaFamily(){
    cout<<"Gamma Family has been created"<<endl;
}
vec GammaFamily::linkfun(vec mu){
    return (1/mu);
}
vec GammaFamily::linkinv(vec eta){
    return (1/eta);
}
vec GammaFamily::varfun(vec mu){
    return(pow(mu,2));
}
vec GammaFamily::dmudeta (vec eta){
    return(-1/pow(eta,2));
}
vec GammaFamily::loglik_vec(vec y, vec eta){
    return this->logdensity(y,eta); 
} 
vec GammaFamily::dloglik(vec y, vec eta){
   vec result(1);
   result<<this->deta(y,eta);
   return result;
}
vec GammaFamily::density(vec y, vec eta){
    return exp(this->logdensity(y,eta));
    //OR
    //return eta%exp(-eta%y);
}

vec GammaFamily::logdensity(vec y, vec eta){
    return log(eta)-eta%y;    
}
double GammaFamily::deta(vec y, vec eta){
    return sum(1/eta-y);
}
