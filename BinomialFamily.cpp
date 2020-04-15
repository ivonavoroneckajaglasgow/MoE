#include "BinomialFamily.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double BinomialFamily::linkfun(double mu){
    return log(mu/(1-mu));
}

