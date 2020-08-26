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
  return 0;
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
  return 0;
}

double ExpertModel::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
  return 0;
}