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
//vec phi(x.size());
//phi.ones();
double logsigma_sq=5.8;
vec mu_beta("0 0"); //prior mean
vec vars("0.01 0.01"); //prior Sigma
mat Sigma_beta=diagmat(vars);

NormalFamily* NF= new NormalFamily();

vec beta0=NF->initialiseBeta(y,X,logsigma_sq);
beta0.print("Intial Beta:");
vec betahat1=NF->findBeta(y,X,logsigma_sq);
betahat1.print("Estimate no prior:");
vec betahat2=NF->findBeta(y,X,logsigma_sq,mu_beta,Sigma_beta);
betahat2.print("Estimate with prior:");
vec betaprop1=NF->proposeBeta(beta0,y,X,logsigma_sq);
betaprop1.print("Proposal no prior:");
vec betaprop2=NF->proposeBeta(beta0,y,X,logsigma_sq,mu_beta,Sigma_beta);
betaprop2.print("Proposal with prior:");



return 0;
}

