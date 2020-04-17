#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "ExpertModel.h"
#include "Family.h"

using namespace std;
using namespace arma;

// g(mu)=X'B=eta
// mu=g_inv(eta)=g_inv(X'B)

class GLMExpert: public ExpertModel{
    public:
    vec y; //observed data
    mat eta; //eta=X'B
    int var; //variance parameter
    Family aFamily; //family object
    vec loglik(vec x); //returns the log-likelihood of the model for a vector of observations x
    vec dloglik(vec x, string param); //returns the derivative of log-likelihood wrt to param of the model for observation x
};