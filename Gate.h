#ifndef MOE_GATE_H
#define MOE_GATE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"
#include "Expert.h"
#include "Node.h"

using namespace std;
using namespace arma;


class Gate: public Node{
    public:
    Gate(); //constructor
    vec gamma; //parameter associated with this gate length rp
    mat Omega; //prior variance for gamma
    void addChild(Node* aChild);
    void printChildren();
    void printDescendants();
    void printTerminalNodes();
    double loglik(mat z, mat pi);
    double loglik(mat X, vec gamma, mat z);
    mat pi_calculator(mat X, vec gamma);
    mat pi_calculator2(mat X, vec gamma);
    mat pi_calculator3(mat X, vec gamma);
    vec score(mat X, mat z, mat pi);
    mat hessian(mat X, mat pi);
    vec findGammaMLEChol(mat X, mat z, mat Omega, mat* L); //Cholesky decomposition
    vec findGammaMLEChol(mat X, mat z, mat Omega);         //Cholesky decomposition
    vec findGammaMLEQR(mat X, mat z, mat Omega, mat* R);   //QR decomposition
    vec findGammaMLEQR(mat X, mat z, mat Omega);           //QR decomposition
    vec findGammaMLE(mat X, mat z, mat Omega);             //no decomposition
    mat makeAchol(vec pi);
    mat makeAchol2 (vec pi);
    mat getXout(mat X, mat pi);
    vec getZeta(mat z, mat pi);
    vec findGammaQR(mat X, mat z, mat Omega, mat* R);
    vec findGammaQR(mat X, mat z, mat Omega);
    vec proposeGamma(vec gammaold, mat X, mat z, mat Omega);
    vec logmvndensity(vec response, vec mean, mat Sigma);
    mat getRowSumsMat(mat A);
};

#endif //MOE_GATE_H