#ifndef MOE_FAMILY_H
#define MOE_FAMILY_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

#include "ExpertModel.h"

//This is a superclass of all Family objects 
//Most functions are virtual and overwritten at subclass levels

class GLMModel:public ExpertModel {
    public:
GLMModel();
virtual vec linkfun(vec mu); //link function for a vector of values of mu
virtual vec linkinv(vec eta); //the inverse of the link function
virtual vec dlinkfun(vec mu);
virtual vec varfun(vec mu); //the variance as function of the mean
virtual vec dmudeta (vec eta);// derivative dmu/deta
virtual vec loglik_vec(vec y, vec eta, double logsigma_sq); //returns the log-likelihood of the model 
virtual vec dloglik(vec y, vec eta, double logsigma_sq); //returns the derivative of log-likelihood 
virtual vec density(vec y, vec eta, double logsigma_sq); // returns density function
virtual vec logdensity(vec y, vec eta, double logsigma_sq); // returns log density function
virtual double deta(vec y, vec eta, double logsigma_sq);// returns the derivative of log likelihood wrt to eta
vec etafun(mat X, vec beta);
virtual double a(double phi); //a(phi) function in the exponential family expression for the family
virtual vec V(vec theta); //b''(theta) function in the exponential family expression for the family
vec initialiseBeta(vec y, mat X, double logsigma_sq);
vec findBeta(vec y, mat X, mat* R, double logsigma_sq);
vec findBeta(vec y, mat X, double logsigma_sq);
vec findBeta(vec y, mat X, mat* R, double logsigma_sq, vec mu_beta, mat Sigma_beta);
vec findBeta(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta);
vec updateBeta(vec betaold, vec y, mat X, double logsigma_sq);
vec updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta);
vec logmvndensity(vec response, vec mean, mat Sigma);
virtual double updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n);
};


#endif //MOE_FAMILY_H