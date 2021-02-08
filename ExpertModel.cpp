#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include <armadillo>

#include "ExpertModel.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Expert Model:: Expert Model object
 * 
 */
ExpertModel::ExpertModel(){
    cout<<"Expert Model has been created."<<endl;
}

/**
 * @brief Log likelihood function value (summed over observations)
 * 
 * @param y observation vector
 * @param eta predictor XBeta
 * @param logsigma_sq variance parameter
 * @return double log-likelihood function value
 */
double ExpertModel::loglik(vec y, vec eta, double logsigma_sq){
 return sum(this->loglik_vec(y,eta,logsigma_sq));
}

/**
 * @brief Log likelihood function
 * Handled at specific Family level
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec vector containing log likelihood function value for each observation
 */
vec ExpertModel::loglik_vec(vec y, vec eta, double logsigma_sq){
return 0;
} 
  
/**
 * @brief Derivatives of log likelihood wrt to eta and sigma^2
 * Handled at specific Family level
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec derivatives of log likelihood wrt to eta and sigma^2
 */
vec ExpertModel::dloglik(vec y, vec eta, double logsigma_sq){
return 0;
} 

/**
 * @brief Density function
 * Handled at specific Family level
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec normal density 
 */
vec ExpertModel::density(vec y, vec eta, double logsigma_sq){
return 0;
} 

/**
 * @brief Density on a log scale
 * Handled at specific Family level
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec density on a log scale
 */
vec ExpertModel::logdensity(vec y, vec eta, double logsigma_sq){
  return 0;
}
    
/**
 * @brief Derivative of log likelihood wrt to eta
 * Handled at specific Family level
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return double derivative of log likelihood wrt to eta
 */
double ExpertModel::deta(vec y, vec eta, double logsigma_sq){
  return 0;
}

/**
 * @brief Calculates predictor Xbeta
 * 
 * @param X design matrix
 * @param beta beta vector
 * @return vec predictor
 */
vec ExpertModel::etafun(mat X, vec beta){
return X*beta;
}

/**
 * @brief Estimates beta 
 * Handled at GLM Model and Normal Model
 * @param y response vector
 * @param X design matrix
 * @param R pointer to R matrix in QR decomposition
 * @param logsigma_sq variance parameter
 * @return vec beta estimate
 */
vec ExpertModel::findBeta(vec y, mat X, mat* R, double logsigma_sq){
  return 0;
}

/**
 * @brief Estimates beta 
 * Handled at GLM Model and Normal Model
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @return vec beta estimate
 */
vec ExpertModel::findBeta(vec y, mat X, double logsigma_sq){
  return 0;
}

/**
 * @brief Estimates beta Bayesian approach
 * Handled at GLM Model and Normal Model
 * @param y response vector
 * @param X design matrix
 * @param R pointer to R matrix in QR decomposition
 * @param logsigma_sq variance parameter
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @return vec beta estimate
 */
vec ExpertModel::findBeta(vec y, mat X, mat* R, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  return 0;
}

/**
 * @brief Estimates beta Bayesian approach
 * Handled at GLM Model and Normal Model
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @return vec beta estimate
 */
vec ExpertModel::findBeta(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  return 0;
}

/**
 * @brief Updates beta
 * Handled at GLM Model and Normal Model
 * @param betaold current beta value
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @return vec update value of beta
 */
vec ExpertModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq){
  return 0;
}

/**
 * @brief Updates beta Bayesian approach
 * Handled at GLM Model and Normal Model
 * @param betaold current beta value
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @return vec updated value of beta
 */
vec ExpertModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  return 0;
}

/**
 * @brief Multivariate normal density on a log scale
 * 
 * @param response response vector
 * @param mean mean vector
 * @param Sigma variance-covariance matrix
 * @return vec multivariate normal density on a log scale
 */
vec ExpertModel::logmvndensity(vec response, vec mean, mat Sigma){
   int k = static_cast<int>(Sigma.n_cols);
   //return 1/(pow(2*M_PI,k/2)*sqrt(det(Sigma)))*exp(-0.5*(response-mean).t()*Sigma.i()*(response-mean)); - not log scale
   return -k/2*log(2*M_PI)-0.5*log(det(Sigma))-0.5*(response-mean).t()*Sigma.i()*(response-mean);
}

/**
 * @brief Multivariate normal density on a log scale
 * 
 * @param response response vector
 * @param mean mean vector
 * @param R pointer to R matrix from the QR decomposition
 * @return vec multivariate normal density on a log scale
 */
vec ExpertModel::logmvndensity(vec response, vec mean, mat* R){
int k=static_cast<int>((*R).n_rows);
//return -k/2*log(2*M_PI)+0.5*sum(log(pow((*R).diag(),2)))-0.5*(response-mean).t()*((*R).t()*(*R))*(response-mean);
mat result;
result=-k/2*log(2*M_PI)+sum(log((*R).diag()))-0.5*sum(pow((*R)*(response-mean),2));
return vectorise(result);
}

/**
 * @brief MLE estimate of log(sigma_sq)
 * Handled at Normal Model
 * @param y response vector
 * @param X design matrix
 * @return double log(sigma_sq) estimate
 */
double ExpertModel::findLogSigmaSq(vec y, mat X){
  return 0;
}

/**
 * @brief Marginal posterior distribution of response y
 * As per Gory details
 * Handled at Normal Model
 * @param y response vector
 * @param X design matrix
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @return double marginal posterior distribution of response y
 */
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
        return 1/randg( distr_param(a,1/b));
    }
    double alpha1=static_cast<double>(a+n/2);
    double alpha2=static_cast<double>(b+sum(pow(y-X*beta,2))/2);
    double sigma_new=1/randg( distr_param(alpha1,1/alpha2)); 
    return log(sigma_new);
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

