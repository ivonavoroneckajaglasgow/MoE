#ifndef MOE_NormalModel_H
#define MOE_NormalModel_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "ExpertModel.h"

using namespace std;
using namespace arma;

//Normal Expert
//f(y)=beta_0+beta_1*x+e
//e~N(0,sigma^2)

class NormalModel: public ExpertModel{
public:
   NormalModel(); //a constructor 
   vec loglik_vec(vec y, vec eta, double logsigma_sq); //returns the log-likelihood 
   vec dloglik(vec y, vec eta, double logsigma_sq);    //returns the vector of derivatives of log-likelihood wrt all parameters
   vec density(vec y, vec eta, double logsigma_sq);    //normal density
   vec logdensity(vec y, vec eta, double logsigma_sq); //log-normal density
   double deta(vec y, vec eta, double logsigma_sq);    //derivative of log-likelihod wrt to eta
   double dsigma (vec y, vec eta, double logsigma_sq); //derivative of the log-likelihood wrt to sigma^2
   vec findBeta(vec y, mat X, double dummy); //MLE estimate for beta
   double findLogSigmaSq(vec y, mat X);      //MLE estimate for sigma
   vec updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta); //Updates beta by drawing from posterior
   double logMarginalPosteriorY(vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b);//Marginal posterior distributio of the response y
   };

#endif //MOE_NormalModel_H