#define _USE_MATH_DEFINES
#define EPS 1e-9

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "GammaFamily.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Gamma Family:: Gamma Family object
 * 
 */
GammaFamily::GammaFamily(){
    cout<<"Gamma Family has been created"<<endl;
}

/**
 * @brief Gamma link function
 * 
 * @param mu mean 
 * @return vec Gamma link
 */
vec GammaFamily::linkfun(vec mu){
    return (1/mu);
}

/**
 * @brief Binomial inverse link function
 * 
 * @param eta predictor
 * @return vec Binomial inverse link
 */
vec GammaFamily::linkinv(vec eta){
    return (1/eta);
}

/**
 * @brief Derivative of the link function wrt mu
 * 
 * @param mu mean
 * @return vec derivative of the link function wrt mu
 */
vec GammaFamily::dlinkfun(vec mu){
return -1/pow(mu,2);
}

/**
 * @brief Variance function
 * 
 * @param mu mean
 * @return vec variance function
 */
vec GammaFamily::varfun(vec mu){
    return(pow(mu,2));
}

/**
 * @brief Derivative of the mean function wrt to the predictor
 * 
 * @param eta predictor
 * @return vec derivative of the mean function wrt to the predictor
 */
vec GammaFamily::dmudeta (vec eta){
    return(-1/pow(eta,2));
}

/**
 * @brief Log likelihood function
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq dummy parameter
 * @return vec vector containing log likelihood function value for each observation
 */
vec GammaFamily::loglik_vec(vec y, vec eta, double logsigma_sq){
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
vec GammaFamily::dloglik(vec y, vec eta, double logsigma_sq){
   vec result(1);
   result<<this->deta(y,eta,logsigma_sq);
   return result;
}

/**
 * @brief Exponential density function
 * 
 * @param y response vector
 * @param eta predictor
 * @param logsigma_sq dummy variable
 * @return vec Exponential density function
 */
vec GammaFamily::density(vec y, vec eta, double logsigma_sq){
    return exp(this->logdensity(y,eta,logsigma_sq));
    //OR
    //return eta%exp(-eta%y);
}

/**
 * @brief Exponential log density function
 * 
 * @param y response vector
 * @param eta predictor
 * @param logsigma_sq dummy variable
 * @return vec Exponential log density function
 */
vec GammaFamily::logdensity(vec y, vec eta, double logsigma_sq){
    return log(eta)-eta%y;    
}

/**
 * @brief Derivative of log likelihood wrt to eta
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq dummy parameter
 * @return double derivative of log likelihood wrt to eta
 */
double GammaFamily::deta(vec y, vec eta, double logsigma_sq){
    return sum(1/eta-y);
}

/**
 * @brief a(phi) function in the exponential family expression for the family
 * 
 * @param phi variable phi in the exponential family expression for the family
 * @return double a(phi)
 */
double GammaFamily::a(double phi){
     return -1;
}

/**
 * @brief b''(theta) function in the exponential family expression for the family
 * 
 * @param theta variable theta in the exponential family expression for the family
 * @return vec b''(theta)
 */
vec GammaFamily::V(vec theta){
    return -1/pow(theta,2);
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
double GammaFamily::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
    return 0;
}