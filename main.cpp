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

vec p_calculator(vec gamma, mat X){
    vec helper=exp(X*gamma);
    return helper/(1+helper);
}

vec etafun(mat X, vec beta){
return X*beta;
}

double transformSigma(double logsigma_sq){
    return exp(logsigma_sq);
}

vec logdensity(vec y, vec eta, double logsigma_sq){ // this should include all experts, now it is just one
   double std;
   double logsigma_tr=transformSigma(logsigma_sq);
   std=sqrt(logsigma_tr);
   return -0.5*(log(2*M_PI)+log(logsigma_tr))-pow(y-eta,2)/(2*logsigma_tr);
}

double loglik(vec y, vec beta, mat X, double logsigma_sq, vec p){
    vec eta=etafun(X,beta);
    double l=sum(log(p)+logdensity(y,eta,logsigma_sq));
    return l;
}

vec logpriordiff_gamma(vec gamma, vec gamma_star, mat gamma_Sigma){
  return -0.5*gamma_star.t()*gamma_Sigma.i()*gamma_star+0.5*gamma.t()*gamma_Sigma.i()*gamma;
}

 vec gammaUpdate(vec gamma, vec y, mat X, vec beta, double logsigma_sq, vec sigma_proposal, mat gamma_Sigma){
    vec p=p_calculator(gamma,X);
    vec v(2,fill::randn);
    vec u=sigma_proposal%v;
    vec gamma_star(gamma.size());
    //gamma_star[0]=gamma[0]+gamma[0]/gamma[1]*u[1]+u[0];
    //gamma_star[1]=gamma[1]+u[1];
    gamma_star[0]=gamma[0]+u[0];
    gamma_star[1]=gamma[1]+u[1];
    vec p_star=p_calculator(gamma_star,X);
    double l=loglik(y,beta,X,logsigma_sq,p);
    double l_star=loglik(y,beta,X,logsigma_sq,p_star);
    double prior_diff=sum(logpriordiff_gamma(gamma,gamma_star,gamma_Sigma));
    double acceptance=l_star-l+prior_diff;
    double m=randu();
    bool   accept=log(m)<acceptance;
    if(accept==1) return gamma_star;
    if(accept==0) return gamma;
}

vec initialiseGamma(vec x){ //only works for a nx2 X
    std::random_device                      rand_dev;
    std::mt19937                            generator(rand_dev());
    std::uniform_real_distribution<double>  distr(x.min(), x.max());
    double x_star=distr(generator);
    vec gamma(2);
    std::random_device rd; 
    std::mt19937       rnd_gen (rd());
    std::exponential_distribution<> rng (0.5);
    double b=rng(rnd_gen);
    gamma[1]=b;
    gamma[0]=-b*x_star;
    return gamma;
}

int main(){

vec y("0.1 12 2.3 4.7 0.5 16 7.3 8 0.9 0.1");
vec x("1 2 3 4 5 6 7 8 9 10");
mat X(10,2);
X.col(0).ones();
X.col(1)=x;
double logsigma_sq=5.8;
vec gamma("4 3");
vec sigma_proposal("0.5 0.5");
vec beta("5 -0.1");
vec gamma_star("1 2");
vec gamma_mu ("1 3");
vec vars("2 2");
mat gamma_Sigma=diagmat(vars);

vec gamma0=initialiseGamma(x);
gamma0.print("Initial gamma:");
vec gamma_new=gammaUpdate(gamma0,y,X,beta,logsigma_sq,sigma_proposal,gamma_Sigma);
gamma_new.print("New gamma:");
return 0;
}

