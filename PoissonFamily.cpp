#define _USE_MATH_DEFINES
#define EPS 1e-9

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
  vec result(eta.size());
  for (int i=0; i<eta.size();i++) 
  result[i]=max(exp(eta[i]),EPS);
  return result;
}

vec PoissonFamily::loglik_vec(vec y, vec eta){
 return this->logdensity(y,eta);
} 

vec PoissonFamily::dloglik(vec y, vec eta){
   vec result(1);
   result<<this->deta(y,eta);
   return result;
} 

vec PoissonFamily::density(vec y, vec eta){
 return exp(this->logdensity(y,eta));
}

vec PoissonFamily::logdensity(vec y, vec eta){
  vec lambda(eta.size()); 
  lambda=exp(eta);
  return -lambda+y%eta-lgamma(y+1);
}

double PoissonFamily::deta(vec y, vec eta){
return sum(y-exp(eta));
}

