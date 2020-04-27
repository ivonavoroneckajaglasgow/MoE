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

}
vec Family::linkinv(vec eta){

}
vec Family::varfun(vec mu){

}
vec Family::dmudeta (vec eta){

}
vec Family::loglik_vec(vec y, vec eta, double logsigma_sq){

}
vec Family::dloglik(vec y, vec eta, double logsigma_sq){

}
vec Family::density(vec y, vec eta, double logsigma_sq){

}
vec Family::logdensity(vec y, vec eta, double logsigma_sq){
    
}
double Family::deta(vec y, vec eta, double logsigma_sq){

}