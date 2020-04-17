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

class NormalExpert: public ExpertModel{
public:
   vec y; //observed data
   vec beta; //a vector of current beta=(intercept,slope) estimates
   double sigma_sq; //a current value of variance sigma^2
   NormalExpert(); //a constructor 
   double loglik(vec x, vec y, vec beta, double sigma_sq); //returns the log-likelihood of the model for a vector of observations x
   double dloglik(vec x, vec y, vec beta, double sigma_sq, string which); //returns the derivative of log-likelihood wrt to param of the model for observation x
   vec dnorm(vec x, vec y, vec beta, double sigma_sq);//normal density
   vec dnorm_log(vec x, vec y, vec beta, double sigma_sq);//wrapper for returning derivatives
   vec dbeta(vec x, vec y, vec beta, double sigma_sq);
   private:
   double transformSigma(double sigma);//Transforms sigma to a log scale
   vec getMu(vec x, vec beta); // Gets a mu from beta and x

};