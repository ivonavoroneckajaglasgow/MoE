#define _USE_MATH_DEFINES
#define EPS 1e-5
#define SigmaMultiple 2

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "GLMModel.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new GLMModel::GLMModel object
 * 
 */
GLMModel::GLMModel(){
 cout<<"GLMModel has been created."<<endl;
}

/**
* @brief Link function
* Handled at specific Family level
* @param mu mean
* @return vec mean
*/
vec GLMModel::linkfun(vec mu){
return 0;
}

/**
 * @brief Inverse link function
 * Handled at specific Family level
 * @param eta predictor
 * @return vec predictor
 */
vec GLMModel::linkinv(vec eta){
return 0;
}

/**
 * @brief Derivative of the link function
 * Handled at specific Family level
 * @param mu mean
 * @return vec derivative of the link function
 */
vec GLMModel::dlinkfun(vec mu){
    return 0;
}

/**
 * @brief Variance function
 * Handled at specific Family level
 * @param mu mean
 * @return vec variance function
 */
vec GLMModel::varfun(vec mu){
return 0;
}

/**
 * @brief Derivative of the mean function wrt to the predictor
 * Handled at specific Family level
 * @param eta predictor
 * @return vec derivative of the mean function wrt to the predictor
 */
vec GLMModel::dmudeta (vec eta){
return 0;
}

/**
 * @brief Log likelihood function
 * Handled at specific Family level
 * @param y response vector
 * @param eta predictor vector
 * @param logsigma_sq variance parameter
 * @return vec vector containing log likelihood function value for each observation
 */
vec GLMModel::loglik_vec(vec y, vec eta, double logsigma_sq){
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
vec GLMModel::dloglik(vec y, vec eta, double logsigma_sq){
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
vec GLMModel::density(vec y, vec eta, double logsigma_sq){
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
vec GLMModel::logdensity(vec y, vec eta, double logsigma_sq){
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
double GLMModel::deta(vec y, vec eta, double logsigma_sq){
return 0;
}

/**
 * @brief a(phi) function in the exponential family expression for the family
 * Handled at specific Family level
 * @param phi variable phi in the exponential family expression for the family
 * @return double a(phi)
 */
double GLMModel::a(double phi){
    return 0;
}

/**
 * @brief b''(theta) function in the exponential family expression for the family
 * Handled at specific Family level
 * @param theta variable theta in the exponential family expression for the family
 * @return vec b''(theta)
 */
vec GLMModel::V(vec theta){
    return 0;
}

/**
 * @brief Estimates beta using the IWLS algorithm
 * 
 * @param y response vector
 * @param X design matrix
 * @param R pointer to R matrix in QR decomposition
 * @param logsigma_sq variance parameter
 * @return vec beta estimate
 */
vec GLMModel::findBeta(vec y, mat X, mat* R, double logsigma_sq){
    vec beta;
    beta=this->initialiseBeta(y,X, logsigma_sq);
    mat Q;
for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z= eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(logsigma_sq)*pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    qr_econ(Q,*R,diagmat(wsqrt)*X);
    vec beta=solve(*R,Q.t()*diagmat(wsqrt)*Z);
    if(all(abs(beta-beta_old)<(EPS,EPS*abs(beta)).max())) break;
}
return beta;
}

/**
 * @brief Estimates beta using the IWLS algorithm
 * Wrapper for the above function
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @return vec beta estimate
 */
vec GLMModel::findBeta(vec y, mat X, double logsigma_sq){
    mat R;
    return this->findBeta(y,X,&R,logsigma_sq);
}

/**
 * @brief Estimates beta using the IWLS algorithm Bayesian approach
 * 
 * @param y response vector
 * @param X design matrix
 * @param R pointer to R matrix in QR decomposition
 * @param logsigma_sq variance parameter
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @return vec estimated beta
 */
vec GLMModel::findBeta(vec y, mat X, mat* R, double logsigma_sq, vec mu_beta, mat Sigma_beta){
    vec beta;
    beta=this->initialiseBeta(y,X,logsigma_sq);
    mat sqrtSigma_beta=sqrtmat_sympd(Sigma_beta);
    mat invSigma_beta=Sigma_beta.i();
    mat sqrt_invSigma_beta=sqrtmat_sympd(invSigma_beta);
    mat Q;
if(X.n_rows==0){
    *R=chol(invSigma_beta);
    return mu_beta;
}
for (int i=0; i<100; i++){
    vec beta_old=beta;
    vec eta=this->etafun(X,beta);
    vec mu=this->linkinv(eta);
    vec Z=eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(logsigma_sq)*pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    mat X_star=join_cols(diagmat(wsqrt)*X,sqrt_invSigma_beta);
    qr_econ(Q,*R,X_star);
    vec Z_star=join_cols(wsqrt%Z,sqrt_invSigma_beta*mu_beta);
    beta=solve(*R,Q.t()*Z_star);
    if(all(abs(beta-beta_old)<(EPS,EPS*abs(beta)).max())) break;
}
return beta;
}

/**
 * @brief Estimates beta using the IWLS algorithm Bayesian approach
 * Wrapper for the function above
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @return vec estimated beta
 */
vec GLMModel::findBeta(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
    mat R;
    return this->findBeta(y,X,&R, logsigma_sq, mu_beta, Sigma_beta);
}

/**
 * @brief Initial step for the IWLS algorithm 
 * 
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @return vec initial value for beta
 */
vec GLMModel::initialiseBeta(vec y, mat X, double logsigma_sq){
    vec mu = y+0.1;
    vec eta= this->linkfun(mu);
    vec Z=eta+(y-mu)%this->dlinkfun(mu);
    vec w=1/(this->a(logsigma_sq)*pow(this->dlinkfun(mu),2)%this->V(mu));
    vec wsqrt=sqrt(w);
    vec beta=solve(diagmat(wsqrt)*X,diagmat(wsqrt)*Z);
    return beta;
}

/**
 * @brief Updates beta
 * 
 * @param betaold current beta value
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @return vec update value of beta
 */
vec GLMModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq){
  mat R;
  vec betahat=this->findBeta(y,X,&R,logsigma_sq);//Uses IWLS to estimate beta
  vec v(betaold.size(),fill::randn);
  vec betanew=betahat+sqrt(SigmaMultiple)*solve(R,v); //take in proposal scale as an argument at some point 
  double density_old=sum(this->logdensity(y,this->etafun(X,betaold),logsigma_sq));
  double density_new=sum(this->logdensity(y,this->etafun(X,betanew),logsigma_sq));
  double proposal_old=sum(this->logmvndensity(betaold,betahat,&R));
  double proposal_new=sum(this->logmvndensity(betanew,betahat,&R));
  double acceptance= density_new-density_old+proposal_old-proposal_new;
  double u=randu();
  bool accept=u<exp(acceptance);
  if(accept==1){
    return betanew;
  }else{
    return betaold; 
  }
} 

/**
 * @brief Updates beta Bayesian approach
 * 
 * @param betaold current beta value
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @return vec updated value of beta
 */
vec GLMModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
  mat R;
  vec betahat=this->findBeta(y,X,&R, logsigma_sq,mu_beta,Sigma_beta);//Uses IWLS to estimate beta
  vec v(betaold.size(),fill::randn);
  vec betanew=betahat+solve(sqrt(SigmaMultiple)*R,v); //take in proposal scale as an argument at some point 
  double density_old=sum(this->logdensity(y,this->etafun(X,betaold),logsigma_sq));
  double density_new=sum(this->logdensity(y,this->etafun(X,betanew),logsigma_sq));
  double proposal_old=sum(this->logmvndensity(betaold,betahat,&R));
  double proposal_new=sum(this->logmvndensity(betanew,betahat,&R));
  double prior_old=sum(this->logmvndensity(betaold,mu_beta,Sigma_beta));
  double prior_new=sum(this->logmvndensity(betanew,mu_beta,Sigma_beta));
  double acceptance=density_new-density_old+proposal_old-proposal_new+prior_new-prior_old;
  double u=randu();
  bool accept=u<exp(acceptance);
  if(accept==1){
       return betanew;
  }else{
    return betaold; 
  }
} 

/**
 * @brief Updates sigma by drawing from the posterior distribution
 * Handled at Expert Model level
 * @param sigma_old old value of variance 
 * @param y respose vector
 * @param X design matrix
 * @param beta beta vector
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param n the number of points in this expert
 * @return double updated value of sigma
 */
double GLMModel::updateSigma(double sigma_old, vec y, mat X, vec beta, double a, double b, int n){
    return 0;
}