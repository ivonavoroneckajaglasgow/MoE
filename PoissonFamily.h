#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Family.h"

using namespace std;
using namespace arma;

//y~Poi(lambda)
//E(y)=mu=lambda
//g(mu)=log(mu)=X'B=eta
//mu=lambda=exp(eta)

class PoissonFamily : public Family {
public:
double linkfun(double mu); //link function for one value of mu
vec linkfun(vec mu); //link function for a vector of values of mu
};