#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Family.h"

using namespace std;
using namespace arma;

//y~Poi(lambda)
//E(y)=mu=lambda
//g(mu)=log(mu)=X'B=eta
//mu=lambda=exp(eta)

class PoissonFamily : public Family {
public:
vec linkfun(vec mu); //link function for a vector of values of mu
vec linkinv(vec eta); //the inverse of the link function
vec varfun(vec mu); //the variance as function of the mean
double dmudeta (vec eta);// derivative dmu/deta
vec loglik_vec(vec y, vec eta, double lambda); //returns the log-likelihood of the model 
vec dloglik(vec y, vec eta, double lambda); //returns the derivative of log-likelihood 
vec density(vec y, vec eta, double lambda); // returns density function
vec logdensity(vec y, vec eta, double lambda); // returns log density function
double deta(vec y, vec eta, double lambda);// returns the derivative of log likelihood wrt to eta
};