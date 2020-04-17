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
    vec x; //explanatory variables
    vec y; //observed data
    mat eta; //eta=X'B
    double var; //variance parameter
    Family aFamily; //family object
    double loglik(vec x); //returns the log-likelihood of the model for a vector of observations x
    double dloglik(vec x, string param); //returns the derivative of log-likelihood wrt to param of the model for observation x
};