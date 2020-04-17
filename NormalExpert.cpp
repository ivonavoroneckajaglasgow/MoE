#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"


#include "NormalExpert.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Normal Expert:: Normal Expert object
 * 
 */
NormalExpert::NormalExpert(){
cout<<"You have created a NormalExpert object."<<endl;
}

/**
 * @brief Mean function mu=X'Beta
 * 
 * @param x vecor of explanatory variables
 * @param beta vector of slope and intercept
 * @return vec mu vector
 */
vec NormalExpert::getMu(vec x, vec beta){
    vec mu;
    mu=beta[0]+beta[1]*x;
    return mu;
}

/**
 * @brief Density of the normal distribution
 * This function automatically transforms sigma by taking an exponential of the input for sigma_sq
 * @param x vector of explanatory variables
 * @param y vector of response variables
 * @param beta beta vector of slope and intercept
 * @param sigma_sq variance parameter 
 * @return vec vector containing density of the normal distribution
 */
vec NormalExpert::dnorm(vec x, vec y, vec beta, double sigma_sq){
   vec mu;
   double sigma_tr;
   sigma_tr=this->transformSigma(sigma_sq);
   mu=this->getMu(x,beta);
   return 1/sqrt(2*M_PI*sigma_tr)*exp(-pow(y-mu,2)/(2*sigma_tr));
}

/**
 * @brief Density of the normal distribution on a log scale
 * This function automatically transforms sigma by taking an exponential of the input for sigma_sq
 * @param x vector of explanatory variables
 * @param y vector of response variables
 * @param beta beta vector of slope and intercept
 * @param sigma_sq variance parameter 
 * @return vec vector containing density of the normal distribution on a log scale
 */
vec NormalExpert::dnorm_log(vec x, vec y, vec beta, double sigma_sq){
   double std;
   double sigma_tr=this->transformSigma(sigma_sq);
   vec mu;
   mu=this->getMu(x,beta);
   std=sqrt(sigma_tr);
   return -0.5*(log(2*M_PI)+log(sigma_tr))-pow(y-mu,2)/(2*sigma_tr);
}

/**
 * @brief Log likelihood function for the normal expert
 * sigma transform takes place on dnorm and dnorm_log level
 * @param x vector of explanatory variables
 * @param y vector of response variables
 * @param beta beta vector of slope and intercept
 * @param sigma_sq variance parameter  
 * @return double value for log likelihood function 
 */
double NormalExpert::loglik(vec x, vec y, vec beta, double sigma_sq){
//return sum(log(this->dnorm(x, beta, sigma_sq)));
return sum(this->dnorm_log(x,y,beta,sigma_sq));
}

/**
 * @brief derivative of the log likelihood function wrt to the intercept or the slope
 * 
 * @param x vector of explanatory variables
 * @param y vector of response variables
 * @param beta vector of slope and intercept
 * @param sigma_sq variance parameter   
 * @param which a choice of intercept (beta0) or slope (beta1)
 * @return value of the derivative of the log likelihood function wrt to the intercept or the slope
 */
double NormalExpert::dloglik(vec x, vec y, vec beta, double sigma_sq, string which){
vec result;
result=this->dbeta(x,y,beta,sigma_sq);
if(which=="beta0") return result[0];
if(which=="beta1") return result[1];
}

/**
 * @brief Transforms supplied value of sigma squared by taking an exponential of it
 * 
 * @param sigma_sq variance parameter
 * @return double transformed value of the variance parameter
 */
double NormalExpert::transformSigma(double sigma){
    return exp(sigma);
}

/**
 * @brief derivative of the log likelihood function wrt to the intercept and the slope
 * 
 * @param x vector of explanatory variables
 * @param y vector of response variables
 * @param beta vector of slope and intercept
 * @param sigma_sq variance parameter   
 * @return vector of the derivative of the log likelihood function wrt to the intercept and the slope
 */
vec NormalExpert::dbeta(vec x, vec y, vec beta, double sigma_sq){
    vec mu;
    vec result;
    double sigma_tr;

    mu=this->getMu(x,beta);
    sigma_tr=this->transformSigma(sigma_sq);

    result<<sum((y-mu)/sigma_tr)<<sum((y-mu)/sigma_tr%x);
    
    return result;
}