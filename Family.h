#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

//This is a superclass of all Family objects 
//Most functions are virtual and overwritten at subclass levels

class Family {
    public:
Family();
virtual vec linkfun(vec mu); //link function for a vector of values of mu
virtual vec linkinv(vec eta); //the inverse of the link function
virtual vec varfun(vec mu); //the variance as function of the mean
virtual vec dmudeta (vec eta);// derivative dmu/deta
virtual vec loglik_vec(vec y, vec eta); //returns the log-likelihood of the model 
virtual vec dloglik(vec y, vec eta); //returns the derivative of log-likelihood 
virtual vec density(vec y, vec eta); // returns density function
virtual vec logdensity(vec y, vec eta); // returns log density function
virtual double deta(vec y, vec eta);// returns the derivative of log likelihood wrt to eta
};


