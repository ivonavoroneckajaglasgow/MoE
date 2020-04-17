#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Family.h"

using namespace std;
using namespace arma;

//y~Bern(1,p)
//E(y)=mu=p
//g(mu)=g(p)=log(p/(1-p))=X'B=eta
//mu=p=exp(eta)/(1+exp(eta))=1/(1+exp(-eta))

class BinomialFamily : public Family {
public:
double linkfun(double mu); //link function for one value of mu
vec linkfun(vec mu); //link function for a vector of values of mu
};