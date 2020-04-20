#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "ExpertModel.h"
#include "Family.h"

using namespace std;
using namespace arma;

// g(mu)=XB=eta
// mu=g_inv(eta)=g_inv(XB)

class GLMExpert: public ExpertModel{
    public:
    Family aFamily; //family object
    virtual vec loglik_vec(vec y, vec eta); //returns the log-likelihood of the model 
    virtual vec dloglik(vec y, vec eta); //returns the derivative of log-likelihood 
    virtual vec density(vec y, vec eta); // returns density function
    virtual vec logdensity(vec y, vec eta); // returns log density function
    virtual double deta(vec y, vec eta);// returns the derivative of log likelihood wrt to eta
};