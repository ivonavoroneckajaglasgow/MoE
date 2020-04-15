#include "PoissonFamily.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double PoissonFamily::linkfun(double mu){
    return log(mu);
}