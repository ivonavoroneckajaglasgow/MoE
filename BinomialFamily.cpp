#define _USE_MATH_DEFINES
#define EPS 1e-9

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "BinomialFamily.h"

using namespace std;
using namespace arma;

BinomialFamily::BinomialFamily(){
    cout<<"Binomial Family has been created."<<endl;
}

vec BinomialFamily::linkfun(vec mu){
  return log(mu/(1-mu));
}

vec BinomialFamily::linkinv(vec eta){
  return exp(eta)/(1+exp(eta)); 
  //OR
  //return 1/(1+exp(-eta)); 
}

vec BinomialFamily::varfun(vec mu){
   return mu%(1-mu); 
}

vec BinomialFamily::dmudeta (vec eta){
    return exp(eta)/(pow(1+exp(eta),2));
}

vec BinomialFamily::loglik_vec(vec y, vec eta){
    return this->logdensity(y,eta);
}

vec BinomialFamily::dloglik(vec y, vec eta){
   vec result(1);
   result<<this->deta(y,eta);
   return result;
}

vec BinomialFamily::density(vec y, vec eta){
return exp(this->logdensity(y,eta));
//OR
/*     vec theta(eta.size());
    theta=this->linkinv(eta);
    vec first(eta.size());
    vec second(eta.size());
    for(int i=0;i<eta.size();i++){
        first[i]=pow(theta[i],y[i]);
        second[i]=pow(1-theta[i],1-y[i]);
    }
    return first%second; */
}

vec BinomialFamily::logdensity(vec y, vec eta){
    vec theta(eta.size());
    theta=this->linkinv(eta);
    return y%log(theta)+(1-y)%log(1-theta);
}

double BinomialFamily::deta(vec y, vec eta){
    return sum(1/(exp(eta)+1));
}