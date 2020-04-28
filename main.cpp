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

int main(){
mat X(3,2);
X="1 3; 1 2; 1 5";
X.print("My matrix: ");
mat Xt;
Xt=X.t();
Xt.print("My matrix transposed: ");
cout<<"Matrix X has "<<X.n_rows<<" rows and "<<X.n_cols<<" columns."<<endl;
cout<<"Matrix Xt has "<<Xt.n_rows<<" rows and "<<Xt.n_cols<<" columns."<<endl;
mat XtX(2,2);
//XtX=Xt * X; //THIS DOES NOT WORK

mat A = randu<mat>(5,10);
mat C = randu<mat>(10,5);
//mat U = A * C; //THIS DOES NOT WORK

}

