#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Family.h"

using namespace std;
using namespace arma;

Family::Family(){
 cout<<"Family has been created."<<endl;
}

vec Family::linkfun(vec mu){
return 0;
}
vec Family::linkinv(vec eta){
return 0;
}
vec Family::varfun(vec mu){
return 0;
}
vec Family::dmudeta (vec eta){
return 0;
}
vec Family::loglik_vec(vec y, vec eta, double logsigma_sq){
return 0;
}
vec Family::dloglik(vec y, vec eta, double logsigma_sq){
return 0;
}
vec Family::density(vec y, vec eta, double logsigma_sq){
return 0;
}
vec Family::logdensity(vec y, vec eta, double logsigma_sq){
    return 0;
}
double Family::deta(vec y, vec eta, double logsigma_sq){
return 0;
}

vec Family::findBeta(vec y, mat X, vec beta){
return 0;
}