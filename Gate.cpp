#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#define _USE_MATH_DEFINES
#define EPS 1e-16

#include "Gate.h"

using namespace std;
using namespace arma;

Gate::Gate(){
    cout<<"Gate has been created."<<endl;
}

double Gate::loglik(mat z, mat pi){
    return sum(vectorise(z%log(pi)))+sum((1-sum(z,1))%log(1-sum(pi,1)));
}

mat Gate::pi_calculator(mat X, vec gamma){
    int p=X.n_cols;
    int r=gamma.size()/p;
    mat gamma2=gamma;
    mat gamma_mat(p,r);

    for(int i=0; i<r;i++){
        mat a=gamma2.rows((i*p),((i+1)*p-1));
        gamma_mat.col(i)=a;
    }

    mat helper=exp(X*gamma_mat);
    mat rowsums(helper.n_rows,helper.n_cols);
        for(int i=0; i<helper.n_cols; i++){
            rowsums.col(i)=sum(helper,1);
        }

    mat final=helper/(1+rowsums);
    mat E(final.n_rows,final.n_cols);
    E.fill(sqrt(EPS));

return arma::min(1-E,arma::max(final,E));
}

vec Gate::score(mat X, mat z, mat pi){
    return vectorise(X.t()*(z-pi));
}

mat Gate::hessian(mat X, mat pi){
 int r=pi.n_cols;
 int p=X.n_cols;

 mat H(r*p,r*p);

 for(int i=0; i<r; i++){
     for(int j=0; j<r; j++){
         mat a;
         if(i==j) a=-pi.col(i)%(1-pi.col(i));
         if(i!=j) a= pi.col(i)%pi.col(j);
        mat a2(X.n_rows,X.n_cols);
            for(int i=0; i<X.n_cols; i++){
            a2.col(i)=a;
        }
    H.submat(span(i*p,(i+1)*p-1),span(j*p,(j+1)*p-1))=X.t()*(a2%X); 
   }
}
 
return H;
}

vec Gate::findGamma(mat X, mat z, mat Omega){
    vec gamma(X.n_cols*z.n_cols);
    gamma.zeros();
    for(int i=0; i<100; i++){
        mat pi=this->pi_calculator(X,gamma);
        vec gamma_old=gamma;
        gamma=gamma-solve(this->hessian(X,pi)-Omega,this->score(X,z,pi)-Omega*gamma);
    if(all(abs(gamma-gamma_old)<(sqrt(EPS),sqrt(EPS)*abs(gamma)).max())) break;
    }
    return gamma;
}

mat Gate::makeAchol(vec pi){
  mat helper(pi.size(),pi.size());
  for(int i=0; i<pi.size();i++){
      for(int j=0; j<pi.size(); j++){
          helper(i,j)=pi[i]*pi[j];
      }
  }
  return chol(diagmat(pi)-helper);
}

mat Gate::getXout(mat X, mat pi){
    int r=pi.n_cols;
    int n=X.n_rows;
    int p=X.n_cols;
    mat Xout(n*r,p*r);

    for(int j=0; j<n; j++){
        mat C=this->makeAchol(vectorise(pi.row(j)));
        mat helper(p*C.n_cols,C.n_cols);
        for(int i=0; i<p; i++){
             vec helper1=vectorise(X.row(j));
             helper.rows(i*r,(i+1)*r-1)=helper1[i]*C;
        }
        helper.reshape(helper.n_cols,helper.n_rows);
        Xout.rows(j*r,(j+1)*r-1)=helper;
    }

    return Xout;
}

vec Gate::getZeta(mat z, mat pi){
    int r=pi.n_cols;
    int n=pi.n_rows;
    vec zeta(n*r);
    for(int i=0; i<n; i++){
        mat C=this->makeAchol(vectorise(pi.row(i)));
        mat helper=z.row(i)-pi.row(i);
        helper.reshape(helper.n_cols,helper.n_rows);
        zeta(span(i*r,(i+1)*r-1))=solve(C.t(),helper);
}
return zeta;
}

vec Gate::findGammaQR(mat X, mat z, mat Omega){
    int r=z.n_cols;
    int n=X.n_rows;
    int p=X.n_cols;
    mat Omega_chol=chol(Omega);
    vec gamma(p*r);
    gamma.zeros();
    for(int i=0; i<100; i++){
        vec gamma_old=gamma;
        mat pi=this->pi_calculator(X,gamma);
        mat Xout=this->getXout(X,pi);
        vec zeta=this->getZeta(z,pi);
        mat LHS=join_cols(Xout,Omega_chol);
        mat RHS(LHS.n_rows,1);
        RHS.zeros();
        RHS.rows(0,zeta.size()-1)=Xout*gamma+zeta;
        mat Q;
        mat R;
        qr(Q,R,LHS);
        gamma=solve(R,Q.t()*RHS);
        if(all(abs(gamma-gamma_old)<(sqrt(EPS),sqrt(EPS)*abs(gamma)).max())) break;
    }
    return gamma;
}