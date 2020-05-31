#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "NormalModel.h"

using namespace std;
using namespace arma;

vec NormalModel::test(vec a){
    return 0;
}

/**
 * @brief Construct a new Normal Expert:: Normal Expert object
 * 
 */
NormalModel::NormalModel(){
cout<<"Normal Model has been created."<<endl;
}

/**
 * @brief function that transforms log(sigma_sq) to sigma_sq by exponentiating it
 * 
 * @param logsigma_sq sigma squared on a log scale 
 * @return double logsigma_sq exponentiated
 */
double NormalModel::transformSigma(double logsigma_sq){
    return exp(logsigma_sq);
}

/**
 * @brief Log density function for a normal distribution
 * Equivalent of dnorm(...,log=TRUE) in R 
 * @param y observed data
 * @param eta vector XB
 * @param logsigma_sq variance parameter
 * @return vec log density function of observed data y for a normal distribution 
 */
vec NormalModel::logdensity(vec y, vec eta, double logsigma_sq){
   double std;
   double logsigma_tr=this->transformSigma(logsigma_sq);
   std=sqrt(logsigma_tr);
   return -0.5*(log(2*M_PI)+log(logsigma_tr))-pow(y-eta,2)/(2*logsigma_tr);
}

/**
 * @brief Density function for a normal distribution
 * Equivalent of dnorm(...) in R 
 * @param y observed data
 * @param eta vector XB
 * @param logsigma_sq variance parameter
 * @return vec density function of observed data y for a normal distribution 
 */
vec NormalModel::density(vec y, vec eta, double logsigma_sq){
    return exp(this->logdensity(y,eta,logsigma_sq));
}

/**
 * @brief log likelihood function for a normal distribution
 * 
 * @param y observed data
 * @param eta vector XB
 * @param logsigma_sq variance parameter
 * @return vec loglik_vec vector of contributions to the log likelihood of each observation y_i
 */
vec NormalModel:: loglik_vec(vec y, vec eta, double logsigma_sq){
    return this->logdensity(y,eta,logsigma_sq);
}

/**
 * @brief derivative of log-likelihood function wrt eta
 * 
 * @param y observed data
 * @param eta vector XB
 * @param logsigma_sq variance parameter
 * @return double derivative of log-likelihood function wrt eta
 */
double NormalModel::deta(vec y, vec eta, double logsigma_sq){
    double logsigma_tr;
    logsigma_tr=this->transformSigma(logsigma_sq);
    return sum((y-eta)/logsigma_tr);
}

/**
 * @brief derivative of log-likelihood function wrt sigma_sq
 * 
 * @param y observed data
 * @param eta vector XB
 * @param logsigma_sq variance parameter
 * @return double derivative of log-likelihood function wrt to sigma_sq
 */
double NormalModel::dsigma(vec y, vec eta, double logsigma_sq){
    double logsigma_tr;
    logsigma_tr=this->transformSigma(logsigma_sq);
    vec density;
    density=this->density(y,eta,logsigma_sq);
    double helper;
    helper=sum(-density*1/(2*logsigma_tr)%(1-pow(y-eta,2)/logsigma_tr));
    double sigma_sq;
    sigma_sq=logsigma_tr*helper;
    return sigma_sq;
}

/**
 * @brief derivative of log-likelihood function wrt to eta and sigma_sq
 * 
 * @param y observed data
 * @param eta vector XB
 * @param logsigma_sq variance parameter
 * @return vec vector containing derivative value of the log-likelihood function wrt eta and sigma_sq (in that order)
 */
vec NormalModel::dloglik(vec y, vec eta, double logsigma_sq){
    vec result;
    result << this->deta(y,eta,logsigma_sq)<<this->dsigma(y,eta,logsigma_sq);
    return result;
}

vec NormalModel::findBeta(vec y, mat X){
mat Xt(X.n_cols,X.n_rows);
Xt=X.t();
mat XtX(X.n_rows,X.n_rows);
XtX=Xt*X;
mat XtXinv(X.n_rows,X.n_rows);
XtXinv=XtX.i();
return XtXinv*Xt*y;
}