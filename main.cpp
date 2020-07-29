#define _USE_MATH_DEFINES
#define EPS 1e-16

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"
#include <chrono>

using namespace std;
using namespace arma;
using namespace std::chrono; 

#include "Gate.h"
#include "NormalFamily.h"

int main(){
cout<<"Structure in C++:"<<endl;

mat X={{0.5,5,87,11},
       {14,1,0.3,0.01},
       {17,7,9.3,11.6}};
int n=X.n_rows;
int p=X.n_cols;
mat z={{0,1},
       {0,0},
       {1,0}};

int r=z.n_cols;

vec diagonals(X.n_cols*z.n_cols);
diagonals.fill(0.001);
mat Omega=diagmat(diagonals);
vec gamma("0.120139084, -1.812376850,  0.151582984, -1.119221005,  0.001908206,  1.188518494, -0.505343855, -0.099234393");

Gate* G1= new Gate();
G1->name="G1";
Gate* G2= new Gate();
G2->name="G2";
Gate* G3= new Gate();
G3->name="G3";
Gate* G4= new Gate();
G4->name="G4";
Expert* E1= new Expert();
E1->name="E1";
Expert* E2= new Expert();
E2->name="E2";
Expert* E3= new Expert();
E3->name="E3";
Expert* E4= new Expert();
E4->name="E4";
Expert* E5= new Expert();
E5->name="E5";

G1->addChild(G2);
G1->addChild(G3);
G2->addChild(E1);
G2->addChild(G4);
G3->addChild(E2);
G3->addChild(E3);
G4->addChild(E4);
G4->addChild(E5);

G1->printChildren();
G4->printChildren();

E4->printParent();
G3->printParent();

cout<<"Descendants of G1:"<<endl;
G1->printDescendants();
cout<<"Terminal nodes of G1:"<<endl;
G1->printTerminalNodes();

cout<<"Descendants of G2:"<<endl;
G2->printDescendants();
cout<<"Terminal nodes of G2:"<<endl;
G2->printTerminalNodes();

G1->gamma=gamma;
G1->Omega=Omega;

vec gammahat=G1->findGammaMLE(X,z,Omega);
gammahat.print("G1 gamma hat:");
vec gammanew=G1->proposeGamma(gamma,X,z,Omega);
gammanew.print("G1 gamma new:");

NormalFamily* NF= new NormalFamily();
// E1->expertmodel=NF; how do I do this?

vec y("1,2,3");
double logsigmasq=1;

vec beta=NF->findBeta(y,X,logsigmasq);
beta.print("betahat:");

mat R;
vec betahat=NF->findBeta(y,X,&R,logsigmasq);//Uses IWLS to estimate beta
betahat.print("betahat:");
vec v(beta.size(),fill::randn);
v.print("v:");
cout<<"R:"<<R.n_rows<<"x"<<R.n_cols<<endl;
R.print("R:");
//Issue occurs here:
//vec betanew=betahat+solve(sqrt(2)*R,v);
//betanew.print("betanew:");


//vec betanew=NF->proposeBeta(beta,y,X,logsigmasq);
//betanew.print("betanew:");

return 0; 
}