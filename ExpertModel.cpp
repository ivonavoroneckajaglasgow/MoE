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

//vec ExpertModel::findBeta(vec y, mat X, vec phi){
 //  return 0; 
//}