#define _USE_MATH_DEFINES
#define EPS 1e-16

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"
#include <chrono>
#include <string> 
#include <algorithm> 
#include <sstream> 
#include <iterator> 

using namespace std;
using namespace arma;
using namespace std::chrono; 

#include "Gate.h"
#include "NormalFamily.h"
#include "Data.h"


int main(){
cout<<"Structure in C++:"<<endl;

Gate* G1= new Gate();
G1->name="G1";
G1->Parent=NULL;
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
Expert* E6= new Expert();
E6->name="E6";

G1->addChild(G2);
G1->addChild(G3);
G1->addChild(E1);
G2->addChild(E2);
G2->addChild(E3);
G3->addChild(E4);
G3->addChild(G4);
G4->addChild(E5);
G4->addChild(E6);
cout<<"EXAMPLE STARTS HERE:"<<endl;
cout<<"First I am going to issue left to right ID down the whole tree."<<endl;
G1->issueIDLR();

int n=10;
int p=2;

mat X(n,p,fill::randn);
vec y(n, fill::randn);

cout<<"Next, create some data with "<< n<<" points and "<<p<<" explanatory variables."<<endl;

vector<Node*> all_experts=G1->getTerminalNodes();
vector<Node*> z_final(n);

cout<<"Now randomly assign the "<<n <<" points to the experts in the tree."<<endl;

for(int i=0;i<n;i++){
   int sub=rand() % all_experts.size();
   z_final[i]=all_experts[sub];
   cout<<"Point "<<i<<" assigned to "<<z_final[i]->name<<endl;
}

cout<<"Now, I will estimate gammas for the gates based on the point allocation:"<<endl;

mat z_G1=G1->getZ(z_final);
mat z_G2=G2->getZ(z_final);
mat z_G3=G3->getZ(z_final);
mat z_G4=G4->getZ(z_final);

mat Omega;

G1->gamma=G1->findGammaMLE(G1->subsetX(X,G1->getPointIndices(z_final)),z_G1,Omega);
cout<<"G1 gamma: "<<G1->gamma<<endl;
G2->gamma=G1->findGammaMLE(G2->subsetX(X,G2->getPointIndices(z_final)),z_G2,Omega);
cout<<"G2 gamma: "<<G2->gamma<<endl;
G3->gamma=G3->findGammaMLE(G3->subsetX(X,G3->getPointIndices(z_final)),z_G3,Omega);
cout<<"G3 gamma: "<<G3->gamma<<endl;
G4->gamma=G4->findGammaMLE(G4->subsetX(X,G4->getPointIndices(z_final)),z_G4,Omega);
cout<<"G4 gamma: "<<G4->gamma<<endl;


cout<<"Next, estimate betas for experts based on the point allocations."<<endl;
cout<<"Set Expert Model to be Normal Family for all experts."<<endl;
cout<<"For simplicity, keep logsigmasq="<<1<<" for all normal experts."<<endl;

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;
E5->expertmodel=NF;
E6->expertmodel=NF;

double logsigma_sq=1;
vec mu_beta(X.n_cols);
mu_beta.zeros();
vec diags(X.n_cols);
diags.fill(0.0001);
mat Sigma_beta=diagmat(diags);

E1->beta=E1->expertmodel->findBeta(E1->subsetY(y,E1->getPointIndices(z_final)),E1->subsetX(X,E1->getPointIndices(z_final)),logsigma_sq,mu_beta,Sigma_beta);
cout<<"Beta for E1: "<<E1->beta<<endl;
E2->beta=E2->expertmodel->findBeta(E2->subsetY(y,E2->getPointIndices(z_final)),E2->subsetX(X,E2->getPointIndices(z_final)),logsigma_sq,mu_beta,Sigma_beta);
cout<<"Beta for E2: "<<E2->beta<<endl;
E3->beta=E3->expertmodel->findBeta(E3->subsetY(y,E3->getPointIndices(z_final)),E3->subsetX(X,E3->getPointIndices(z_final)),logsigma_sq,mu_beta,Sigma_beta);
cout<<"Beta for E3: "<<E3->beta<<endl;
E4->beta=E4->expertmodel->findBeta(E4->subsetY(y,E4->getPointIndices(z_final)),E4->subsetX(X,E4->getPointIndices(z_final)),logsigma_sq,mu_beta,Sigma_beta);
cout<<"Beta for E4: "<<E4->beta<<endl;
E5->beta=E5->expertmodel->findBeta(E5->subsetY(y,E5->getPointIndices(z_final)),E5->subsetX(X,E5->getPointIndices(z_final)),logsigma_sq,mu_beta,Sigma_beta);
cout<<"Beta for E5: "<<E5->beta<<endl;
E6->beta=E6->expertmodel->findBeta(E6->subsetY(y,E6->getPointIndices(z_final)),E6->subsetX(X,E6->getPointIndices(z_final)),logsigma_sq,mu_beta,Sigma_beta);
cout<<"Beta for E6: "<<E6->beta<<endl;

E1->logsigma_sq=1;
E2->logsigma_sq=1;
E3->logsigma_sq=1;
E4->logsigma_sq=1;
E5->logsigma_sq=1;
E6->logsigma_sq=1;

cout<<"Now we can update chosen betas and gammas."<<endl;

E1->beta=E1->expertmodel->updateBeta(E1->beta,E1->subsetY(y,E1->getPointIndices(z_final)),E1->subsetX(X,E1->getPointIndices(z_final)),logsigma_sq,mu_beta,Sigma_beta);
cout<<"New E1 beta:"<<E1->beta<<endl;
G1->gamma=G1->updateGamma(G1->gamma,X,z_G1,Omega);
cout<<"New G1 gamma:"<<G1->gamma<<endl;

cout<<"Finally, update z_final:"<<endl;
cout<<"Choosing to update from G1 downwards,but could do any other gate"<<endl;

vector<Node*> z_new=G1->updateZ(y,X,z_final);

cout<<"Look at old allocations vs new:"<<endl;

for(int i=0;i<z_final.size();i++){
       cout<<"Point "<<i<<" before was in "<<z_final[i]->name<<" and now is in "<<z_new[i]->name<<endl;
}

cout<<"Now one can use getZ() to get new matrix z for any gate based on z_new"<<endl;


mat X_new(n,p,fill::randn);
vec yhat=G1->predict(X_new);
yhat.print("yhat:");


double y2=5;
double a=1;
double b=2;

double sigma_new=NF->updateSigma(E1->logsigma_sq,y,X,E1->beta,a,b,y.size());
cout<<"Sigma update test: "<<sigma_new<<endl;


vector<Node*> z_MCMC=G1->MCMC(10,y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_final);

for(int i=0;i<z_final.size();i++){
  cout<<"Before: "<<z_final[i]->name<<". After: "<<z_MCMC[i]->name<<endl;
}


//string E1JSON=E1->createJSON();
//cout<<E1JSON<<endl;
//string E4JSON=E4->createJSON();
//cout<<E4JSON<<endl;
//string G1JSON=G1->createJSON();
//cout<<G1JSON<<endl;
//cout<<"STRING TEST"<<endl;

//string test=G1->createJSON3();
//cout<<"test: "<<test<<endl;

// for(int i=0;i<test.size();i++){
//   cout<<"Character "<<i<<": "<<test[i]<<endl;
// }

//vec test1=G1->logmvndensity(E1->beta,mu_beta,Sigma_beta);
//test1.print("test 1: ");

return 0;
}