//
// Created by 2022175V on 15/08/2019.
//

#ifndef MOE_3_GATEPARAMETERS_H
#define MOE_3_GATEPARAMETERS_H

#include "armadillo"
using namespace arma;

struct GateParameters {
    vec gamma;
    vec prior_gamma_mean;
    mat prior_gamma_var;
};


#endif //MOE_3_GATEPARAMETERS_H
