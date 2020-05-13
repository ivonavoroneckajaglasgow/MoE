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
 vec x("1 2 3 4 5 6 7 8 9 10");
 vec y("12 2.3 3.4 0.4 25 36 17 8.5 9.1 0.10");
 vec w("0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1");
 vec wsqrt=sqrt(w);
 vec phi(x.size());
 phi.ones();
 mat X(10,2);
 X.col(0).ones();
 X.col(1)=x;
 NormalFamily* NF=new NormalFamily();
 vec betahat=NF->findBeta(y,X,phi);
 betahat.print("Estimated betahat:");
 mat Q,R;
 qr_econ(Q,R,diagmat(wsqrt)*X);
 R.print("My R matrix from QR:");
 vec betaproposal=NF->proposeBeta(betahat,R);
 betaproposal.print("Proposed beta:");
}

