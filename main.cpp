#define _USE_MATH_DEFINES
#define EPS 1e-9

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

#include "Expert.h"
#include "ExpertModel.h"
#include "NormalModel.h"
#include "GLMModel.h"
#include "Family.h"
#include "NormalFamily.h"
#include "PoissonFamily.h"
#include "BinomialFamily.h" 
#include "GammaFamily.h"

int main(){
 vec x("1 2 3 4 5 6 7 8 9 10");
 //x.print("vector x:");
 mat X(10,2);
 X.col(0).ones();
 X.col(1)=x;
 //X.print("Matrix X:");
 vec w("0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0");
 mat W;
 W=diagmat(w);
 //W.print("Matrix W:");
vec z("10 4 1 8 2 9 6 3 7 5");
vec betahat=solve(X.t()*W*X,X.t()*W*z);
//betahat.print("betahat:");
//Choleski decomposition:
mat C=chol(X.t()*W*X);
//C.print("C:");
mat result(10000,2);
result.ones();
for (int i=0;i<result.n_rows;i++){
  vec v(2,fill::randn);
  vec a=betahat+solve(C,v);
  rowvec rowa(2);
  rowa<<a[0]<<a[1];
  result.row(i)=rowa;
}
cout<<"Cholesky Decomposition:"<<endl;
cout<<"Column Means:"<<endl;
cout<< mean(result, 0);
rowvec betahatrow(2);
betahatrow<<betahat[0]<<betahat[1];
betahatrow.print("Compare to betahat:");

cout<<"Estimated Sigma:"<<endl;
cout<<cov(result);
cout<<"Compare to:"<<endl;
cout<<(X.t()*W*X).i()<<endl;

//QR Decomposition
cout<<"QR Decomposition:"<<endl;
mat Q,R;
mat Wsqrt;
Wsqrt=sqrt(W);
qr_econ(Q,R,Wsqrt*X);
vec betahat2;
betahat2=solve(R,Q.t()*Wsqrt*z);
mat result2(10000,2);
result2.ones();
for (int i=0;i<result2.n_rows;i++){
  vec v2(2,fill::randn);
  vec a2=betahat+solve(R,v2);
  rowvec rowa2(2);
  rowa2<<a2[0]<<a2[1];
  result2.row(i)=rowa2;
}

cout<<"Column Means:"<<endl;
cout<< mean(result2, 0);
rowvec betahatrow2(2);
betahatrow2<<betahat2[0]<<betahat2[1];
betahatrow2.print("Compare to betahat:");

cout<<"Estimated Sigma:"<<endl;
cout<<cov(result2);
cout<<"Compare to:"<<endl;
cout<<(X.t()*W*X).i()<<endl;

}

