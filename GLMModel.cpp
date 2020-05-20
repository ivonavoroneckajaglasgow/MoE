#define _USE_MATH_DEFINES
#define EPS 1e-5
#define SigmaMultiple 2

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

vec GLMModel::findBeta(vec y, mat X, vec phi, mat* R){
    vec beta;
    beta=this->initialiseBeta(y,X,phi);
    mat Q;
for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z= eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(phi)%pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    qr_econ(Q,*R,diagmat(wsqrt)*X);
    vec beta=solve(*R,Q.t()*diagmat(wsqrt)*Z);
    if(all(abs(beta-beta_old)<(EPS,EPS*abs(beta)).max())) break;
}
return beta;
}

vec GLMModel::findBeta(vec y, mat X, vec phi){
    mat R;
    return this->findBeta(y,X,phi,&R);
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
  vec newbeta=beta+solve(sqrt(SigmaMultiple)*R,v);
  //need to add accept reject step based on likelihood
  return newbeta; 
}

