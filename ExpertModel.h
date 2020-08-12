#ifndef MOE_EXPERTMODEL_H
#define MOE_EXPERTMODEL_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

//ExpertModel is a superclass and doesn't do anything itself
//Most functions are virtual and overwritten at subclasses level

class ExpertModel{
    public:
    ExpertModel(); //constructor
    double loglik(vec y, vec eta, double logsigma_sq); //log likelihood function (returns one value summed up over observations)
    virtual vec loglik_vec(vec y, vec eta, double logsigma_sq); //returns the log-likelihood as a vector of likelihood contribution of each observation
    virtual vec dloglik(vec y, vec eta, double logsigma_sq); //returns the derivative of log-likelihood wrt to all parameters
    virtual vec density(vec y, vec eta, double logsigma_sq); //returns the density function
    virtual vec logdensity(vec y, vec eta, double logsigma_sq); // returns the log density function
    virtual double deta(vec y, vec eta, double logsigma_sq);// returns the derivative of the log likelihood wrt eta
    //virtual vec findBeta(vec y, mat X, vec phi);
    virtual vec etafun(mat X, vec beta);
    virtual vec initialiseBeta(vec y, mat X, double logsigma_sq);
    virtual vec findBeta(vec y, mat X, mat* R, double logsigma_sq);
    virtual vec findBeta(vec y, mat X, double logsigma_sq);
    virtual vec findBeta(vec y, mat X, mat* R, double logsigma_sq, vec mu_beta, mat Sigma_beta);
    virtual vec findBeta(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta);
    virtual vec proposeBeta(vec betaold, vec y, mat X, double logsigma_sq);
    virtual vec proposeBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta);
    virtual vec logmvndensity(vec response, vec mean, mat Sigma);
};

#endif //MOE_EXPERTMODEL_H