#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "ExpertModel.h"

using namespace std;
using namespace arma;

ExpertModel::ExpertModel(){
    cout<<"Expert Model has been creaeted."<<endl;
}

double ExpertModel::loglik(vec y, vec eta){

}

vec ExpertModel::loglik_vec(vec y, vec eta){

} 
  
vec ExpertModel::dloglik(vec y, vec eta){

} 

vec ExpertModel::density(vec y, vec eta){

} 

vec ExpertModel::logdensity(vec y, vec eta){

}
    
double ExpertModel::deta(vec y, vec eta){

}