#define _USE_MATH_DEFINES
#define EPS 1e-16

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

#include "Gate.h"

int main(){
cout<<"SETTING UP GATE"<<endl;

mat X={{0.1,0.4,0.7,1},
       {0.2,0.5,0.8,1.1},
       {0.3,0.6,0.9,1.2}};
X.print("Design matrix X:");
int n=X.n_rows;
cout<<"Number of observations n="<<n<<endl;
int p=X.n_cols;
cout<<"Number of covariates p="<<p<<endl;

mat z={{0,1},
       {0,0},
       {1,0}};

int r=z.n_cols;
cout<<"Number of splits (excluding ref class) r="<<r<<endl;
z.print("Matrix of allocations z:");
cout<<"Each row corresponds to an observation and each column to a split."<<endl;

Gate* G= new Gate();

vec gamma("-0.84085548,1.38435934,-1.25549186,0.07014277,1.71144087,-0.60290798,-0.47216639,-0.63537131");
gamma.print("Gating parameters gamma:");
cout<<"Length of gamma (rp): "<<gamma.size()<<endl;

mat pi=G->pi_calculator(X,gamma);
pi.print("pi:");

double loglik=G->loglik(z,pi);
cout<<"Log-likelihood:"<<loglik<<endl;

vec score=G->score(X,z,pi);
score.print("Score:");

mat H=G->hessian(X,pi);
H.print("Hessian:");

vec diagonals(X.n_cols*z.n_cols);
diagonals.fill(0.001);
mat Omega=diagmat(diagonals);
Omega.print("Omega (inverse of a prior varcov matrix):");

vec gamma_hat=G->findGamma(X,z,Omega);
gamma_hat.print("gammahat:");

// mat test=G->makeAchol(vectorise(pi.row(1)));
// test.print("test:");

mat Xout=G->getXout(X,pi);
cout<<"Xout is a matrix of size nr x pr: "<<Xout.n_rows<<" x "<<Xout.n_cols<<endl;
Xout.print("Xout:");

//perform a couple random checks and compare to R
cout<<Xout(4,7)<<endl;
cout<<Xout(3,0)<<endl;
cout<<Xout(2,6)<<endl;
cout<<Xout(3,3)<<endl;

vec zeta=G->getZeta(z,pi);
zeta.print("zeta:");

vec gamma_hatQR= G->findGammaQR(X,z,Omega);
gamma_hatQR.print("gammahat QR:");

cout<<"Differences in estimates with and without QR:"<<endl;
cout<<abs(gamma_hat-gamma_hatQR)<<endl;

vec gammanew=G->proposeGamma(gamma,X,z,Omega);
gammanew.print("One MCMC step update:");

return 0; 
}