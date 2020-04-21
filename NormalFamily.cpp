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

vec NormalFamily::varfun(vec mu){
    vec varvec(mu.size());
    varvec.ones();
    return varvec; 
}

vec NormalFamily::dmudeta(vec eta){
    //L said this has to be one number, but in R it is a vector - double-check this
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




