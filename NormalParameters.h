//
// Created by 2022175V on 15/08/2019.
//

#ifndef MOE_3_NORMALPARAMETERS_H
#define MOE_3_NORMALPARAMETERS_H

#include "armadillo"
using namespace arma;

struct NormalParameters{
    vec beta;
    double sigma;
    vec prior_beta;
    mat V;
    double A;
    double B;
};

#endif //MOE_3_NORMALPARAMETERS_H
