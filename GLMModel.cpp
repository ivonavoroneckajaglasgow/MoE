#define _USE_MATH_DEFINES
#define EPS 1e-5

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "GLMModel.h"

using namespace std;
using namespace arma;

GLMModel::GLMModel(){
 cout<<"GLMModel has been created."<<endl;
}

vec GLMModel::linkfun(vec mu){
return 0;
}
vec GLMModel::linkinv(vec eta){
return 0;
}
vec GLMModel::dlinkfun(vec mu){
    return 0;
}
vec GLMModel::varfun(vec mu){
return 0;
}
vec GLMModel::dmudeta (vec eta){
return 0;
}
vec GLMModel::loglik_vec(vec y, vec eta, double logsigma_sq){
return 0;
}
vec GLMModel::dloglik(vec y, vec eta, double logsigma_sq){
return 0;
}
vec GLMModel::density(vec y, vec eta, double logsigma_sq){
return 0;
}
vec GLMModel::logdensity(vec y, vec eta, double logsigma_sq){
    return 0;
}
double GLMModel::deta(vec y, vec eta, double logsigma_sq){
return 0;
}

vec GLMModel::etafun(mat X, vec beta){
return X*beta;
}

vec GLMModel::a(vec phi){
    return 0;
}

vec GLMModel::V(vec theta){
    return 0;
}

vec GLMModel::findBeta(vec y, mat X, vec phi){
    vec beta;
    beta=this->initialiseBeta(y,X,phi);

for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z= eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(phi)%pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    mat Q,R;
    qr_econ(Q,R,diagmat(wsqrt)*X);
    vec beta=solve(R,Q.t()*diagmat(wsqrt)*Z);
    //vec beta=solve(diagmat(wsqrt)*X,diagmat(wsqrt)*Z); //line without QR decomp
    if(all(abs(beta-beta_old)<(EPS,EPS*abs(beta)).max())) break;
}
return beta;
}

mat GLMModel::findR(vec y, mat X, vec phi){
    vec beta;
    beta=this->initialiseBeta(y,X,phi);

for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z= eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(phi)%pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    mat Q,R;
    qr_econ(Q,R,diagmat(wsqrt)*X);
    vec beta=solve(R,Q.t()*diagmat(wsqrt)*Z);
    //vec beta=solve(diagmat(wsqrt)*X,diagmat(wsqrt)*Z); //line without QR decomp
    if(all(abs(beta-beta_old)<(EPS,EPS*abs(beta)).max())) break;
}
return 0;
//return R; will not work  
}

vec GLMModel::initialiseBeta(vec y, mat X, vec phi){
    vec mu = y+0.1;
    vec eta= this->linkfun(mu);
    vec Z=eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(phi)%pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    vec beta=solve(diagmat(wsqrt)*X,diagmat(wsqrt)*Z);
    return beta;
}

vec GLMModel::proposeBeta(vec beta, mat R){
  vec v(2,fill::randn);
  vec newbeta=beta+solve(R,v);
  //need to add accept reject step based on likelihood
  return newbeta; 
}

