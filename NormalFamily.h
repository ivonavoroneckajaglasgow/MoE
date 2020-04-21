#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Family.h"

using namespace std;
using namespace arma;

//y~N(mu,sigma^2)
//E(y)=mu
//g(mu)=mu=X'B=eta
//mu=eta

class NormalFamily : public Family {
public:
NormalFamily();
vec linkfun(vec mu); //link function for a vector of values of mu
vec linkinv(vec eta); //the inverse of the link function
vec varfun(vec mu); //the variance as function of the mean
vec dmudeta (vec eta);// derivative dmu/deta
vec loglik_vec(vec y, vec eta, double var); //returns the log-likelihood of the model 
vec dloglik(vec y, vec eta, double var); //returns the derivative of log-likelihood 
vec density(vec y, vec eta, double var); // returns density function
vec logdensity(vec y, vec eta, double var); // returns log density function
double deta(vec y, vec eta, double var);// returns the derivative of log likelihood wrt to eta
double dsigma (vec y, vec eta, double logsigma_sq); //derivative of the log-likelihood wrt to sigma^2
private:   
double transformSigma(double logsigma_sq);//Transforms sigma to a log scale
};