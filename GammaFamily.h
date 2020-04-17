#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Family.h"

using namespace std;
using namespace arma;

//y~Gamma(alpha,beta)
//E(y)=mu=alpha/beta
//g(mu)=-1/mu=X'B=eta
//mu=-1/eta

class GammaFamily : public Family {
public:
double linkfun(double mu); //link function for one value of mu
vec linkfun(vec mu); //link function for a vector of values of mu
};