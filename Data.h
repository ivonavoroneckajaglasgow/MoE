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
    mat X;
    vec y;
    Data();
};

#endif //MOE_DATA_H