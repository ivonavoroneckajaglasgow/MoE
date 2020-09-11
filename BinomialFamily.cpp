#define _USE_MATH_DEFINES
#define EPS 1e-9

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "BinomialFamily.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Binomial Family:: Binomial Family object
 * 
 */
BinomialFamily::BinomialFamily(){
    cout<<"Binomial Family has been created."<<endl;
}

/**
 * @brief Binomial link function
 * 
 * @param mu mean 
 * @return vec Binomial link
 */
vec BinomialFamily::linkfun(vec mu){
  return log(mu/(1+mu));
}

/**
 * @brief Binomial inverse link function
 * 
 * @param eta predictor
 * @return vec Binomial inverse link
 */
vec BinomialFamily::linkinv(vec eta){
  return exp(eta)/(1+exp(eta)); 
  //OR
  //return 1/(1+exp(-eta)); 
}

/**
 * @brief Derivative of the link function wrt mu
 * 
 * @param mu mean
 * @return vec derivative of the link function wrt mu
 */
vec BinomialFamily::dlinkfun(vec mu){
  return 1/(mu%(1-mu));
}

/**
 * @brief Variance function
 * 
 * @param mu mean
 * @return vec variance function
 */
vec BinomialFamily::varfun(vec mu){
   return mu%(1-mu); 
}

/**
 * @brief Derivative of the mean function wrt to the predictor
 * 
 * @param eta predictor
 * @return vec derivative of the mean function wrt to the predictor
 */
vec BinomialFamily::dmudeta (vec eta){
    return exp(eta)/(pow(1+exp(eta),2));
}

/**
 * @brief Log likelihood function
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq dummy parameter
 * @return vec vector containing log likelihood function value for each observation
 */
vec BinomialFamily::loglik_vec(vec y, vec eta, double logsigma_sq){
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
vec BinomialFamily::dloglik(vec y, vec eta, double logsigma_sq){
   vec result(1);
   result<<this->deta(y,eta,logsigma_sq);
   return result;
}

/**
 * @brief Binomial density function
 * 
 * @param y response vector
 * @param eta predictor
 * @param logsigma_sq dummy variable
 * @return vec Binomial density function
 */
vec BinomialFamily::density(vec y, vec eta, double logsigma_sq){
return exp(this->logdensity(y,eta,logsigma_sq));
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

/**
 * @brief Binomial log density function
 * 
 * @param y response vector
 * @param eta predictor
 * @param logsigma_sq dummy variable
 * @return vec Binomial log density function
 */
vec BinomialFamily::logdensity(vec y, vec eta, double logsigma_sq){
    vec theta(eta.size());
    theta=this->linkinv(eta);
    return y%log(theta)+(1-y)%log(1-theta);
}

/**
 * @brief Derivative of log likelihood wrt to eta
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq dummy parameter
 * @return double derivative of log likelihood wrt to eta
 */
double BinomialFamily::deta(vec y, vec eta, double logsigma_sq){
    return sum(1/(exp(eta)+1));
}

/**
 * @brief a(phi) function in the exponential family expression for the family
 * 
 * @param phi variable phi in the exponential family expression for the family
 * @return double a(phi)
 */
double BinomialFamily::a(double phi){
    return 1;
}

/**
 * @brief b''(theta) function in the exponential family expression for the family
 * 
 * @param theta variable theta in the exponential family expression for the family
 * @return vec b''(theta)
 */
vec BinomialFamily::V(vec theta){
    return exp(theta)/pow(1+exp(theta),2);
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
double BinomialFamily::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
    return 0;
}