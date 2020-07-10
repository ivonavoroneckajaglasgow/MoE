#ifndef MOE_GATE_H
#define MOE_GATE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;


class Gate{
    public:
    Gate(); //constructor
    double loglik(mat z, mat pi);
    mat pi_calculator(mat X, vec gamma);
    vec score(mat X, mat z, mat pi);
    mat hessian(mat X, mat pi);
    vec findGamma(mat X, mat z, mat Omega);
    mat makeAchol(vec pi);
    mat getXout(mat X, mat pi);
    vec getZeta(mat z, mat pi);
    vec findGammaQR(mat X, mat z, mat Omega);
    vec proposeGamma(vec gammaold, mat X, mat z, mat Omega);
    vec logmvndensity(vec response, vec mean, mat Sigma);
};

#endif //MOE_GATE_H