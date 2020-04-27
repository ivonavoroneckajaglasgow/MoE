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
 cout<<"Creating expert:"<<endl;
 Expert* E= new Expert(); //creates expert model as well
 cout<<"Creating normal model:"<<endl;
 ExpertModel* EM=new ExpertModel();
 cout<<"Creating normal expert:"<<endl;
 NormalModel* NE= new NormalModel(); //creates expert model as well
 cout<<"Creating GLM model:"<<endl;
 GLMModel* GE= new GLMModel(); //creates expert model as well
 cout<<"Creating a family:"<<endl;
 Family* F=new Family();
 cout<<"Creating a normal family:"<<endl;
 NormalFamily* NF=new NormalFamily(); //creates a family as well
 cout<<"Creating a poisson family:"<<endl;
 PoissonFamily* PF=new PoissonFamily(); //creates a family as well
cout<<"Creating a binomial family:"<<endl;
 BinomialFamily* BF=new BinomialFamily(); //creates a family as well
}

