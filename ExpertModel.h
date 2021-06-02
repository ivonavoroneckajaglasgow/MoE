#ifndef MOE_EXPERTMODEL_H
#define MOE_EXPERTMODEL_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Node.h"

using namespace std;
using namespace arma;

//ExpertModel is a superclass and doesn't do anything itself
//Most functions are virtual and overwritten at subclasses level

class ExpertModel{
    public:
    ExpertModel(); //constructor
    double loglik(vec y, vec eta, double logsigma_sq);          //log likelihood function (returns one value summed up over observations)
    virtual vec loglik_vec(vec y, vec eta, double logsigma_sq); //returns the log-likelihood as a vector of likelihood contribution of each observation
    virtual vec dloglik(vec y, vec eta, double logsigma_sq);    //returns the derivative of log-likelihood wrt to all parameters
    virtual vec density(vec y, vec eta, double logsigma_sq);    //returns the density function
    virtual vec logdensity(vec y, vec eta, double logsigma_sq); // returns the log density function
    virtual double deta(vec y, vec eta, double logsigma_sq);    // returns the derivative of the log likelihood wrt eta
    vec etafun(mat X, vec beta); //XBeta calculator
    virtual vec findBeta(vec y, mat X, mat* R, double logsigma_sq); //estimates beta + pointer to R in QR
    virtual vec findBeta(vec y, mat X, double logsigma_sq);         //estimates beta
    virtual vec findBeta(vec y, mat X, mat* R, double logsigma_sq, vec mu_beta, mat Sigma_beta);//estimates beta + pointer to R in QR Bayesian
    virtual vec findBeta(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta);        //estimates beta Bayesian
    virtual double findLogSigmaSq(vec y, mat X); //estimates sigma sq on a log scale
    virtual vec updateBeta(vec betaold, vec y, mat X, double logsigma_sq); //updates beta
    virtual vec updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta); //updates beta Bayesian
    vec logmvndensity(vec response, vec mean, mat Sigma); //multivariate normal density with variance-covariance matrix
    vec logmvndensity(vec response, vec mean, mat* R);    //multivariate normal density with R from QR decomposition
    virtual double logMarginalPosteriorY(vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b); //marginal posterior distribution of response y
    double updateSigma(vec y, mat X, vec beta, double a, double b, int n);//updates sigma
    double IG_log(double y, double a, double b);//inverse gamma density on a log scale   
    double transformSigma(double logsigma_sq);  //Transforms sigma to a log scale
    virtual vec findBetaMLE(vec y, mat X);
    virtual double findLogSigmaSqMLE(vec y, mat X, vec betahat);    
    double qSigma(vec y, mat X, vec beta, double logsigma_sq, double a, double b);
    virtual double qBeta(vec y, mat X, vec beta, double logsigma_sq, vec mu_beta, mat Sigma_beta);
};

#endif //MOE_EXPERTMODEL_H