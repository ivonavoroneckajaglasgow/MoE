#ifndef MOE_DATA_H
#define MOE_DATA_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

class Data {
    public:
    mat X; //design matrix
    vec y; //response vector
    Data();//constructor
};

#endif //MOE_DATA_H