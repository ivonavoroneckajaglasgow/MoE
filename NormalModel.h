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
   vec dloglik(vec y, vec eta, double logsigma_sq); //returns the vector of derivatives of log-likelihood wrt all parameters
   vec density(vec y, vec eta, double logsigma_sq);//normal density
   vec logdensity(vec y, vec eta, double logsigma_sq);//log-normal density
   double deta(vec y, vec eta, double logsigma_sq); //derivative of log-likelihod wrt to eta
   double dsigma (vec y, vec eta, double logsigma_sq); //derivative of the log-likelihood wrt to sigma^2
   vec findBeta(vec y, mat X, vec beta);
   private:   
   double transformSigma(double logsigma_sq);//Transforms sigma to a log scale
   };

#endif //MOE_NormalModel_H