#define _USE_MATH_DEFINES
#define EPS 1e-9

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

#include "Expert.h"
#include "ExpertModel.h"
#include "NormalModel.h"
#include "GLMModel.h"
#include "NormalFamily.h"
#include "PoissonFamily.h"
#include "BinomialFamily.h" 
#include "GammaFamily.h"

int main(){
vec y("0.1 12 2.3 4.7 0.5 16 7.3 8 0.9 0.1");
vec x("1 2 3 4 5 6 7 8 9 10");
mat X(10,2);
X.col(0).ones();
X.col(1)=x;
vec w("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1");
vec wsqrt=sqrt(w);
vec phi(x.size());
phi.ones();

NormalFamily* GLM = new NormalFamily();
vec beta=GLM->findBeta(y,X,phi);
beta.print("beta:");
//How do I access R?

mat R;
vec beta2=GLM->findBeta(y,X,phi,&R);
R.print("R:");

}

