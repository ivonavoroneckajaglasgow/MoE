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
    vec x; //explanatory variables
    vec y; //observed data
    ExpertModel(); //constructor
    virtual double loglik(vec x); //returns the log-likelihood of the model for a vector of observations x
    virtual double dloglik(vec x, string param); //returns the derivative of log-likelihood wrt to param of the model for observation x
};