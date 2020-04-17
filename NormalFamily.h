#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Family.h"

using namespace std;
using namespace arma;

//y~N(mu,sigma^2)
//E(y)=mu
//g(mu)=mu=X'B=eta
//mu=eta

class NormalFamily : public Family {
public:
double linkfun(double mu); //link function for one value of mu
vec linkfun(vec mu); //link function for a vector of values of mu
};