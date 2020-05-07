#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "GLMModel.h"

using namespace std;
using namespace arma;

GLMModel::GLMModel(){
 cout<<"GLM Model has been created"<<endl;
} 

vec GLMModel::loglik_vec(vec y, vec eta, double logsigma_sq){
return 0;
}

vec GLMModel::dloglik(vec y, vec eta, double logsigma_sq){
return 0;
}
    
vec GLMModel::density(vec y, vec eta, double logsigma_sq){
return 0;
}
    
vec GLMModel::logdensity(vec y, vec eta, double logsigma_sq){
    return 0;
}
    
double GLMModel::deta(vec y, vec eta, double logsigma_sq){
return 0;
}

vec GLMModel::findBeta(vec y, mat X){
return 0;
}