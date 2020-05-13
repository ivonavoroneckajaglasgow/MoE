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
};

#endif //MOE_EXPERTMODEL_H