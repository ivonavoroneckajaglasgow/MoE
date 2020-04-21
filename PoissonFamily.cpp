#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "PoissonFamily.h"

using namespace std;
using namespace arma;

PoissonFamily::PoissonFamily(){
cout<<"Poisson Family has been created"<<endl;
}

vec PoissonFamily::linkfun(vec mu){
    return log(mu);
}

vec PoissonFamily::linkinv(vec eta){
    return exp(eta);
}

vec PoissonFamily::varfun(vec mu){
    return mu;
}

vec PoissonFamily::dmudeta(vec eta){
    //return exp(eta).max();
}

vec PoissonFamily::loglik_vec(vec y, vec eta, double lambda){

} 
vec PoissonFamily::dloglik(vec y, vec eta, double lambda){

} 
vec PoissonFamily::density(vec y, vec eta, double lambda){

}
vec PoissonFamily::logdensity(vec y, vec eta, double lambda){

}
double PoissonFamily::deta(vec y, vec eta, double lambda){
    
}

