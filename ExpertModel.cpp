#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "ExpertModel.h"

using namespace std;
using namespace arma;

ExpertModel::ExpertModel(){
    cout<<"Expert Model has been created."<<endl;
}

double ExpertModel::loglik(vec y, vec eta, double logsigma_sq){
 return sum(this->loglik_vec(y,eta,logsigma_sq));
}

vec ExpertModel::loglik_vec(vec y, vec eta, double logsigma_sq){
return 0;
} 
  
vec ExpertModel::dloglik(vec y, vec eta, double logsigma_sq){
return 0;
} 

vec ExpertModel::density(vec y, vec eta, double logsigma_sq){
return 0;
} 

vec ExpertModel::logdensity(vec y, vec eta, double logsigma_sq){
  return 0;
}
    
double ExpertModel::deta(vec y, vec eta, double logsigma_sq){
  return 0;
}

vec ExpertModel::etafun(mat X, vec beta){
return X*beta;
}

vec ExpertModel::initialiseBeta(vec y, mat X, double logsigma_sq){
  return 0;
}
vec ExpertModel::findBeta(vec y, mat X, mat* R, double logsigma_sq){
  return 0;
}
vec ExpertModel::findBeta(vec y, mat X, double logsigma_sq){
  return 0;
}
vec ExpertModel::findBeta(vec y, mat X, mat* R, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  return 0;
}
vec ExpertModel::findBeta(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  return 0;
}
vec ExpertModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq){
  return 0;
}
vec ExpertModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  return 0;
}

vec ExpertModel::logmvndensity(vec response, vec mean, mat Sigma){
   int k = Sigma.n_cols;
   //return 1/(pow(2*M_PI,k/2)*sqrt(det(Sigma)))*exp(-0.5*(response-mean).t()*Sigma.i()*(response-mean)); - not log scale
   return -k/2*log(2*M_PI)-0.5*log(det(Sigma))-0.5*(response-mean).t()*Sigma.i()*(response-mean);
}

vec ExpertModel::logmvndensity(vec response, vec mean, mat* R){
int k=(*R).n_rows;
//return -k/2*log(2*M_PI)+0.5*sum(log(pow((*R).diag(),2)))-0.5*(response-mean).t()*((*R).t()*(*R))*(response-mean);
double result=-k/2*log(2*M_PI)+sum(log((*R).diag()))-0.5*sum(pow((*R)*(response-mean),2));
return vectorise(result);
}

double ExpertModel::findLogSigmaSq(vec y, mat X){
  return 0;
}

double ExpertModel::logMarginalPosteriorY(vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b){
  return 0;
}

/**
 * @brief Updates sigma by drawing from the posterior distribution
 * If an expert is empty draws from the prior
 * @param sigma_old old value of variance 
 * @param y respose vector
 * @param X design matrix
 * @param beta beta vector
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param n the number of points in this expert
 * @return double updated value of sigma
 */
double ExpertModel::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
    if(n==0){
        return 1/randg( distr_param(a,b));
    }
    double alpha1=static_cast<double>(a+n/2);
    double alpha2=static_cast<double>(b+sum(pow(y-X*beta,2))/2);
    double sigma_new=1/randg( distr_param(alpha1,alpha2)); 
    return sigma_new;
    // double density_old=sum(this->logdensity(y,this->etafun(X,beta),sigma_old));
    // double density_new=sum(this->logdensity(y,this->etafun(X,beta),sigma_new));
    // double prior_old=sum(this->IG_log(sigma_old,a,b));
    // double prior_new=sum(this->IG_log(sigma_new,a,b));
    // double acceptance=density_new-density_old+prior_new-prior_old;
    // double u=randu();
    // bool accept=u<exp(acceptance);
    // if(accept==1) return sigma_new;
    // if(accept==0) return sigma_old; 
}

/**
 * @brief Inverse Gamma density on a log scale
 * 
 * @param y response vector
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @return double Inverse Gamma density on a log scale
 */
double ExpertModel::IG_log(double y, double a, double b){
    return  a*log(b)-lgamma(a)-(a+1)*log(y)-b/y;
}

/**
 * @brief Function that transforms log(sigma_sq) to sigma_sq by exponentiating it
 * 
 * @param logsigma_sq sigma squared on a log scale 
 * @return double logsigma_sq exponentiated
 */
double ExpertModel::transformSigma(double logsigma_sq){
    return exp(logsigma_sq);
}