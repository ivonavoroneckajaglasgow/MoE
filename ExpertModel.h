#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Expert.h"

using namespace std;
using namespace arma;

//ExpertModel is a superclass and doesn't do anything itself
//Most functions are virtual and overwritten at subclasses level

class ExpertModel: public Expert{
    public:
    ExpertModel(); //constructor
    virtual vec loglik_vec(vec y, vec eta); //returns the log-likelihood as a vector of likelihood contribution of each observation
    virtual vec dloglik(vec y, vec eta); //returns the derivative of log-likelihood wrt to all parameters
    virtual vec density(vec y, vec eta); //returns the density function
    virtual vec logdensity(vec y, vec eta); // returns the log density function
    virtual double deta(vec y, vec eta);// returns the derivative of the log likelihood wrt eta
};