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

}

vec GLMModel::dloglik(vec y, vec eta, double logsigma_sq){

}
    
vec GLMModel::density(vec y, vec eta, double logsigma_sq){

}
    
vec GLMModel::logdensity(vec y, vec eta, double logsigma_sq){
    
}
    
double GLMModel::deta(vec y, vec eta, double logsigma_sq){

}