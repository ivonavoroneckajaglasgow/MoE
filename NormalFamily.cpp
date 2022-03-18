#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "NormalFamily.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Normal Family:: Normal Family object
 * 
 */
NormalFamily::NormalFamily(){
//cout<<"Normal Family has been created"<<endl;
}

/**
 * @brief Identity link function
 * 
 * @param mu mean
 * @return vec mean
 */
vec NormalFamily::linkfun(vec mu){
    return mu;
}

/**
 * @brief Inverse link function
 * 
 * @param eta predictor
 * @return vec predictor
 */
vec NormalFamily::linkinv(vec eta){
    return eta;
}

/**
 * @brief Derivative of the link function
 * 
 * @param mu mean
 * @return vec derivative of the link function
 */
vec NormalFamily::dlinkfun(vec mu){
    vec result(mu.size());
    result.ones();
    return result;
}

/**
 * @brief Variance function
 * 
 * @param mu mean
 * @return vec variance function
 */
vec NormalFamily::varfun(vec mu){
    vec varvec(mu.size());
    varvec.ones();
    return varvec; 
}

/**
 * @brief Derivative of the mean function wrt to the predictor
 * 
 * @param eta predictor
 * @return vec derivative of the mean function wrt to the predictor
 */
vec NormalFamily::dmudeta(vec eta){
    vec dvec(eta.size());
    dvec.ones();
    return dvec;
}

/**
 * @brief Normal density on a log scale
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec normal density on a log scale
 */
vec NormalFamily::logdensity(vec y, vec eta, double logsigma_sq){
   double std;
   //cout<<"eta check:"<<eta[0]<<endl;
   //cout<<"logsigma"<<logsigma_sq<<endl;
   //double logsigma_tr=this->transformSigma(logsigma_sq);
   double logsigma_tr=exp(logsigma_sq);
   std=sqrt(logsigma_tr);
   return -0.5*(log(2*M_PI)+log(logsigma_tr))-pow(y-eta,2)/(2*logsigma_tr);
}

/**
 * @brief Normal density
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec normal density 
 */
vec NormalFamily::density(vec y, vec eta, double logsigma_sq){
    return exp(this->logdensity(y,eta,logsigma_sq));
}

/**
 * @brief Log likelihood function
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec vector containing log likelihood function value for each observation
 */
vec NormalFamily::loglik_vec(vec y, vec eta, double logsigma_sq){
    return this->logdensity(y,eta,logsigma_sq);
}

/**
 * @brief Derivative of log likelihood wrt to eta
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return double derivative of log likelihood wrt to eta
 */
double NormalFamily::deta(vec y, vec eta, double logsigma_sq){
    double logsigma_tr;
    logsigma_tr=this->transformSigma(logsigma_sq);
    
    return sum((y-eta)/logsigma_tr);
}

/**
 * @brief Derivative of the log-likelihood wrt to sigma^2
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return double derivative of log likelihood wrt to sigma^2
 */
double NormalFamily::dsigma(vec y, vec eta, double logsigma_sq){
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
 * @brief Derivatives of log likelihood wrt to eta and sigma^2
 * 
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec derivatives of log likelihood wrt to eta and sigma^2
 */
vec NormalFamily::dloglik(vec y, vec eta, double logsigma_sq){
    vec result;
    result << this->deta(y,eta,logsigma_sq)<<this->dsigma(y,eta,logsigma_sq);
    return result;
}

/**
 * @brief a(phi) function in the exponential family expression for the family
 * 
 * @param phi variable phi in the exponential family expression for the family
 * @return double a(phi)
 */
double NormalFamily::a(double phi){
    //vec result(phi.size());
    //result.ones();
    //return result;
    return phi;
}

/**
 * @brief b''(theta) function in the exponential family expression for the family
 * 
 * @param theta variable theta in the exponential family expression for the family
 * @return vec b''(theta)
 */
vec NormalFamily::V(vec theta){
    vec result(theta.size());
    result.ones();
    return result;
}

double NormalFamily::findLogSigmaSqMLE(vec y, mat X, vec betahat){
    double n=static_cast<int>(y.size());
    //return as_scalar(log(sqrt(1/(n-X.n_cols)*(y-X*betahat).t()*(y-X*betahat))));
    return as_scalar(log(1/(n-X.n_cols)*(y-X*betahat).t()*(y-X*betahat)));
}


vec NormalFamily::findBetaMLE(vec y, mat X){
mat Xt(X.n_cols,X.n_rows);
Xt=X.t();
mat XtX(X.n_rows,X.n_rows);
XtX=Xt*X;
mat XtXinv(X.n_rows,X.n_rows);
XtXinv=XtX.i();
return XtXinv*Xt*y;
}


