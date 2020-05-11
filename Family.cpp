#define _USE_MATH_DEFINES
#define EPS 1e-5

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
vec Family::dlinkfun(vec mu){
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

vec Family::etafun(mat X, vec beta){
return X*beta;
}

vec Family::a(vec phi){
    return 0;
}

vec Family::V(vec theta){
    return 0;
}

vec Family::findBeta(vec y, mat X, vec phi){
    vec beta;
    beta=this->initialiseBeta(y,X,phi);

for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z= eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(phi)%pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    vec beta=solve(diagmat(wsqrt)*X,diagmat(wsqrt)*Z);
    if(any(abs(beta-beta_old)<EPS)) break;
}
return beta;
}

vec Family::initialiseBeta(vec y, mat X, vec phi){
    vec mu = y+0.1;
    vec eta= this->linkfun(mu);
    vec Z=eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(phi)%pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    vec beta=solve(diagmat(wsqrt)*X,diagmat(wsqrt)*Z);
    return beta;
}