
//
// Created by 2022175V on 15/08/2019.
//

#ifndef MOE_3_NORMALEXPERT_H
#define MOE_3_NORMALEXPERT_H

#include "Expert.h"
#include "NormalParameters.h"
#include "armadillo"


using namespace std;
using namespace arma;


class NormalExpert: public Expert {
public:
    NormalParameters parameters;

    NormalExpert(string aName, NormalParameters aParameters);
};


#endif //MOE_3_NORMALEXPERT_H