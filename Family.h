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
virtual double linkfun(double mu); //link function for one value of mu
virtual vec linkfun(vec mu); //link function for a vector of values of mu
virtual vec linkinv(vec mu); //the inverse of the link function
virtual vec varfun(vec mu); //the variance as function of the mean
};
