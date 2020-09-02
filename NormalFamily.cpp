#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "NormalFamily.h"

using namespace std;
using namespace arma;

NormalFamily::NormalFamily(){
cout<<"Normal Family has been created"<<endl;
}

double NormalFamily::transformSigma(double logsigma_sq){
    return exp(logsigma_sq);
}


vec NormalFamily::linkfun(vec mu){
    return mu;
}

vec NormalFamily::linkinv(vec eta){
    return eta;
}

vec NormalFamily::dlinkfun(vec mu){
    vec result(mu.size());
    result.ones();
    return result;
}

vec NormalFamily::varfun(vec mu){
    vec varvec(mu.size());
    varvec.ones();
    return varvec; 
}

vec NormalFamily::dmudeta(vec eta){
    vec dvec(eta.size());
    dvec.ones();
    return dvec;
}

vec NormalFamily::logdensity(vec y, vec eta, double logsigma_sq){
   double std;
   double logsigma_tr=this->transformSigma(logsigma_sq);
   std=sqrt(logsigma_tr);
   return -0.5*(log(2*M_PI)+log(logsigma_tr))-pow(y-eta,2)/(2*logsigma_tr);
}

vec NormalFamily::density(vec y, vec eta, double logsigma_sq){
    return exp(this->logdensity(y,eta,logsigma_sq));
}

vec NormalFamily::loglik_vec(vec y, vec eta, double logsigma_sq){
    return this->logdensity(y,eta,logsigma_sq);
}

double NormalFamily::deta(vec y, vec eta, double logsigma_sq){
    double logsigma_tr;
    logsigma_tr=this->transformSigma(logsigma_sq);
    
    return sum((y-eta)/logsigma_tr);
}

double NormalFamily::dsigma(vec y, vec eta, double logsigma_sq){
    double logsigma_tr;
    logsigma_tr=this->transformSigma(logsigma_sq);
    vec density;
    density=this->density(y,eta,logsigma_sq);
    double helper;
    helper=sum(-density*1/(2*logsigma_tr)%(1-pow(y-eta,2)/logsigma_tr));
    double sigma_sq;
    sigma_sq=logsigma_tr*helper;
    return sigma_sq;
}

vec NormalFamily::dloglik(vec y, vec eta, double logsigma_sq){
    vec result;
    result << this->deta(y,eta,logsigma_sq)<<this->dsigma(y,eta,logsigma_sq);
    return result;
}

double NormalFamily::a(double phi){
    //vec result(phi.size());
    //result.ones();
    //return result;
    return phi;
}

vec NormalFamily::V(vec theta){
    vec result(theta.size());
    result.ones();
    return result;
}

double NormalFamily::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
    if(n==0){
        return 1/randg( distr_param(a,b));
    }
    double alpha1=a+n/2;
    double alpha2=b+sum(pow(y-X*beta,2))/2;
    double sigma_new=1/randg( distr_param(alpha1,alpha2)); 
    double density_old=sum(this->logdensity(y,this->etafun(X,beta),sigma_old));
    double density_new=sum(this->logdensity(y,this->etafun(X,beta),sigma_new));
    double prior_old=sum(this->IG_log(sigma_old,a,b));
    double prior_new=sum(this->IG_log(sigma_new,a,b));
    double acceptance=density_new-density_old+prior_new-prior_old;
    double u=randu();
    bool accept=u<exp(acceptance);
    if(accept==1) return sigma_new;
    if(accept==0) return sigma_old; 
}

double NormalFamily::IG_log(double y, double a, double b){
    return  a*log(b)-lgamma(a)-(a+1)*log(y)-b/y;
}

