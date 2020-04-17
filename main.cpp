#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

#include  "NormalExpert.h"

int main(){

NormalExpert* E = new NormalExpert();

vec beta(2);
double sigma_sq=2;
beta<<2<<1;
vec x(3);
x<<1<<2<<5;
vec y(3);
y<<3<<0.5<<1;

vec test_dnorm;
test_dnorm=E->dnorm(x,y,beta,sigma_sq);
test_dnorm.print("Test dnorm:");

vec test_dnormlog;
test_dnormlog=E->dnorm_log(x,y,beta,sigma_sq);
test_dnormlog.print("Test dnorm_log:");

vec test_loglik;
test_loglik=E->loglik(x,y,beta,sigma_sq);
test_loglik.print("Test loglik:");

vec test_dbeta;
test_dbeta=E->dbeta(x,y,beta,sigma_sq);
test_dbeta.print("Test dbeta:");

double test_beta0;
test_beta0=E->dloglik(x,y,beta,sigma_sq,"beta0");
cout<<"DBeta0:"<<test_beta0<<endl;

double test_beta1;
test_beta1=E->dloglik(x,y,beta,sigma_sq,"beta1");
cout<<"DBeta1:"<<test_beta1<<endl;

return 0;
}

