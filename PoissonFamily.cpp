#define _USE_MATH_DEFINES
#define EPS 1e-9

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "PoissonFamily.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Poisson Family:: Poisson Family object
 * 
 */
PoissonFamily::PoissonFamily(){
cout<<"Poisson Family has been created"<<endl;
}

/**
 * @brief Log link function
 * 
 * @param mu mean 
 * @return vec log link
 */
vec PoissonFamily::linkfun(vec mu){
    return log(mu);
}

/**
 * @brief Log inverse link function
 * 
 * @param eta predictor
 * @return vec log inverse link
 */
vec PoissonFamily::linkinv(vec eta){
    return exp(eta);
}

/**
 * @brief Derivative of the link function wrt mu
 * 
 * @param mu mean
 * @return vec derivative of the link function wrt mu
 */
vec PoissonFamily::dlinkfun(vec mu){
  return 1/mu;
}

/**
 * @brief Variance function
 * 
 * @param mu mean
 * @return vec variance function
 */
vec PoissonFamily::varfun(vec mu){
    return mu;
}

/**
 * @brief Derivative of the mean function wrt to the predictor
 * 
 * @param eta predictor
 * @return vec derivative of the mean function wrt to the predictor
 */
vec PoissonFamily::dmudeta(vec eta){
  vec result(eta.size());
  for (int i=0; i<eta.size();i++) 
  result[i]=max(exp(eta[i]),EPS);
  return result;
}

/**
 * @brief Log likelihood function
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq dummy parameter
 * @return vec vector containing log likelihood function value for each observation
 */
vec PoissonFamily::loglik_vec(vec y, vec eta, double logsigma_sq){
 return this->logdensity(y,eta,logsigma_sq);
} 

/**
 * @brief Derivative of log likelihood wrt to eta
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq dummy parameter
 * @return vec derivative of log likelihood wrt to eta 
 */
vec PoissonFamily::dloglik(vec y, vec eta, double logsigma_sq){
   vec result(1);
   result<<this->deta(y,eta,logsigma_sq);
   return result;
} 

/**
 * @brief Poisson density function
 * 
 * @param y response vector
 * @param eta predictor
 * @param logsigma_sq dummy variable
 * @return vec Poisson density function
 */
vec PoissonFamily::density(vec y, vec eta, double logsigma_sq){
 return exp(this->logdensity(y,eta,logsigma_sq));
}

/**
 * @brief Poisson log density function
 * 
 * @param y response vector
 * @param eta predictor
 * @param logsigma_sq dummy variable
 * @return vec Poisson log density function
 */
vec PoissonFamily::logdensity(vec y, vec eta, double logsigma_sq){
  vec lambda(eta.size()); 
  lambda=exp(eta);
  return -lambda+y%eta-lgamma(y+1);
}

/**
 * @brief Derivative of log likelihood wrt to eta
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq dummy parameter
 * @return double derivative of log likelihood wrt to eta
 */
double PoissonFamily::deta(vec y, vec eta, double logsigma_sq){
return sum(y-exp(eta));
}

/**
 * @brief a(phi) function in the exponential family expression for the family
 * 
 * @param phi variable phi in the exponential family expression for the family
 * @return double a(phi)
 */
double PoissonFamily::a(double phi){
  return 1;
}

/**
 * @brief b''(theta) function in the exponential family expression for the family
 * 
 * @param theta variable theta in the exponential family expression for the family
 * @return vec b''(theta)
 */
vec PoissonFamily::V(vec theta){
  return exp(theta);
}

/**
 * @brief Dummy function that does nothing
 * 
 * @param sigma_old old value of variance 
 * @param y respose vector
 * @param X design matrix
 * @param beta beta vector
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param n the number of points in this expert
 * @return double updated value of sigma
 */
double PoissonFamily::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
    return 0;
}
