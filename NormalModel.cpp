#define _USE_MATH_DEFINES
#define EPS 1e-5
#define SigmaMultiple 2

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "NormalModel.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Normal Expert:: Normal Expert object
 * 
 */
NormalModel::NormalModel(){
cout<<"Normal Model has been created."<<endl;
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

/**
 * @brief MLE estimate for beta
 * 
 * @param y response vector
 * @param X design matrix
 * @param dummy dummy variable taht ensures we can call findBeta from Expert Model
 * @return vec beta estimate
 */
vec NormalModel::findBeta(vec y, mat X, double dummy){
mat Xt(X.n_cols,X.n_rows);
Xt=X.t();
mat XtX(X.n_rows,X.n_rows);
XtX=Xt*X;
mat XtXinv(X.n_rows,X.n_rows);
XtXinv=XtX.i();
return XtXinv*Xt*y;
}

/**
 * @brief MLE estimate of log(sigma_sq)
 * 
 * @param y response vector
 * @param X design matrix
 * @return double log(sigma_sq) estimate
 */
double NormalModel::findLogSigmaSq(vec y, mat X){
    vec betahat=this->findBeta(y,X,1);
    double n=static_cast<int>(y.size());
    return as_scalar(log(1/n*(y-X*betahat).t()*(y-X*betahat)));
}

/**
 * @brief Updates beta by drawing from posterior
 * If an expert is empty draws from the prior
 * @param betaold old value of beta
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance parameter
 * @param mu_beta prior mean of beta
 * @param Sigma_beta prior variance-covariance matrix of beta
 * @return vec new value of beta
 */
 vec NormalModel::updateBeta(vec betaold, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
     mat var;
     vec mean;
     if(X.n_rows==0){
         mean=mu_beta;
         var=Sigma_beta;
     }else{
         mat invSigma_beta=Sigma_beta.i();
         var=(invSigma_beta+X.t()*X).i();
         mean=var*(invSigma_beta*mu_beta+X.t()*y);
    }
     mat R=chol(var);
     vec v(mean.size(),fill::randn);
     vec betanew=mean+sqrt(SigmaMultiple)*solve(R,v);
     return betanew;
    //  double density_old=sum(this->logdensity(y,this->etafun(X,betaold),logsigma_sq));
    //  double density_new=sum(this->logdensity(y,this->etafun(X,betanew),logsigma_sq));
    //  double proposal_old=sum(this->logmvndensity(betaold,mean,&R));
    //  double proposal_new=sum(this->logmvndensity(betanew,mean,&R));
    //  double prior_old=sum(this->logmvndensity(betaold,mu_beta,Sigma_beta));
    //  double prior_new=sum(this->logmvndensity(betanew,mu_beta,Sigma_beta));
    //  double acceptance=density_new-density_old+proposal_old-proposal_new+prior_new-prior_old;
    //  double u=randu();
    //  bool accept=u<exp(acceptance);
    //  if(accept==1) return betanew;
    //  if(accept==0) return betaold;
 }

/**
 * @brief Marginal posterior distribution of response y
 * As per Gory details
 * @param y response vector
 * @param X design matrix
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @return double marginal posterior distribution of response y
 */
double NormalModel::logMarginalPosteriorY(vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b){
    double n=static_cast<double>(X.n_rows);
    mat Sigma_beta_inv=Sigma_beta.i();
    mat Sigma_beta_post=(Sigma_beta_inv+X.t()*X).i();
    vec mu_beta_post=Sigma_beta_post*(Sigma_beta_inv*mu_beta+X.t()*y);
    double sign;
    double detval;
    log_det(detval,sign,Sigma_beta_post);
    double detval2;
    double sign2;
    log_det(detval2,sign2,Sigma_beta);
    double C=a*log(b)+lgamma(a+n/2)+0.5*detval-(n/2)*log(2*M_PI)-lgamma(a)-0.5*detval2;
    double helper=as_scalar(b+0.5*(mu_beta.t()*Sigma_beta_inv*mu_beta+y.t()*y-mu_beta_post.t()*Sigma_beta_post.i()*mu_beta_post));
    //cout<<C-(a+n/2)*log(helper)<<endl;
    return C-(a+n/2)*log(helper);
}