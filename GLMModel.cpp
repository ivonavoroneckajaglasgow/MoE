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

double GLMModel::a(double phi){
    return 0;
}

vec GLMModel::V(vec theta){
    return 0;
}

vec GLMModel::findBeta(vec y, mat X, mat* R, double logsigma_sq){
    vec beta;
    beta=this->initialiseBeta(y,X, logsigma_sq);
    mat Q;
for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z= eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(logsigma_sq)*pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    qr_econ(Q,*R,diagmat(wsqrt)*X);
    vec beta=solve(*R,Q.t()*diagmat(wsqrt)*Z);
    if(all(abs(beta-beta_old)<(EPS,EPS*abs(beta)).max())) break;
}
return beta;
}

vec GLMModel::findBeta(vec y, mat X, double logsigma_sq){
    mat R;
    return this->findBeta(y,X,&R,logsigma_sq);
}

vec GLMModel::findBeta(vec y, mat X, mat* R, double logsigma_sq, vec mu_beta, mat Sigma_beta){
    vec beta;
    beta=this->initialiseBeta(y,X,logsigma_sq);
    mat sqrtSigma_beta=sqrtmat_sympd(Sigma_beta);
    mat invSigma_beta=Sigma_beta.i();
    mat sqrt_invSigma_beta=sqrtmat_sympd(invSigma_beta);
    mat Q;
if(X.n_rows==0){
    *R=chol(invSigma_beta);
    return mu_beta;
}
for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z=eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(logsigma_sq)*pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    mat X_star=join_cols(diagmat(wsqrt)*X,sqrt_invSigma_beta);
    qr_econ(Q,*R,X_star);
    vec Z_star=join_cols(wsqrt%Z,sqrt_invSigma_beta*mu_beta);
    beta=solve(*R,Q.t()*Z_star);
    if(all(abs(beta-beta_old)<(EPS,EPS*abs(beta)).max())) break;
}
return beta;
}

vec GLMModel::findBeta(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
    mat R;
    return this->findBeta(y,X,&R, logsigma_sq, mu_beta, Sigma_beta);
}

vec GLMModel::initialiseBeta(vec y, mat X, double logsigma_sq){
    vec mu = y+0.1;
    vec eta= this->linkfun(mu);
    vec Z=eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(logsigma_sq)*pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    vec beta=solve(diagmat(wsqrt)*X,diagmat(wsqrt)*Z);
    return beta;
}

vec GLMModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq){
  mat R;
  vec betahat=this->findBeta(y,X,&R,logsigma_sq);//Uses IWLS to estimate beta
  vec v(betaold.size(),fill::randn);
  vec betanew=betahat+solve(sqrt(SigmaMultiple)*R,v); //take in proposal scale as an argument at some point 
  double density_old=sum(this->logdensity(y,this->etafun(X,betaold),logsigma_sq));
  double density_new=sum(this->logdensity(y,this->etafun(X,betanew),logsigma_sq));
  double proposal_old=sum(this->logmvndensity(betaold,betahat,&R));
  double proposal_new=sum(this->logmvndensity(betanew,betahat,&R));
  double acceptance= density_new-density_old+proposal_old-proposal_new;
  double u=randu();
  bool accept=u<exp(acceptance);
  if(accept==1) return betanew;
  if(accept==0) return betaold; 
} 

vec GLMModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  mat R;
  vec betahat=this->findBeta(y,X,&R, logsigma_sq,mu_beta,Sigma_beta);//Uses IWLS to estimate beta
  vec v(betaold.size(),fill::randn);
  vec betanew=betahat+solve(sqrt(SigmaMultiple)*R,v); //take in proposal scale as an argument at some point 
  double density_old=sum(this->logdensity(y,this->etafun(X,betaold),logsigma_sq));
  double density_new=sum(this->logdensity(y,this->etafun(X,betanew),logsigma_sq));
  double proposal_old=sum(this->logmvndensity(betaold,betahat,&R));
  double proposal_new=sum(this->logmvndensity(betanew,betahat,&R));
  double prior_old=sum(this->logmvndensity(betaold,mu_beta,Sigma_beta));
  double prior_new=sum(this->logmvndensity(betanew,mu_beta,Sigma_beta));
  double acceptance=density_new-density_old+proposal_old-proposal_new+prior_new-prior_old;
  double u=randu();
  bool accept=u<exp(acceptance);
  if(accept==1) return betanew;
  if(accept==0) return betaold; 
} 

vec GLMModel::logmvndensity(vec response, vec mean, mat Sigma){
   int k = Sigma.n_cols;
   //return 1/(pow(2*M_PI,k/2)*sqrt(det(Sigma)))*exp(-0.5*(response-mean).t()*Sigma.i()*(response-mean)); - not log scale
   return -k/2*log(2*M_PI)-0.5*log(det(Sigma))-0.5*(response-mean).t()*Sigma.i()*(response-mean);
}

vec GLMModel::logmvndensity(vec response, vec mean, mat* R){
int k=(*R).n_rows;
return -k/2*log(2*M_PI)+0.5*sum(log(pow((*R).diag(),2)))-0.5*(response-mean).t()*((*R).t()*(*R))*(response-mean);
}

double GLMModel::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
    return 0;
}