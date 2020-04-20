#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

//This is a superclass of NormalExert and ExpertModel objects 
//Most functions are virtual and overwritten at subclass levels

class Expert{
    public:
    vec y;       //data observed
    vec eta;     // XB
    double var;  // some variance parameter   
    virtual double loglik(vec y, vec eta); //log likelihood function (returns one value summed up over observations)
};