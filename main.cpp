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
#include "NormalModel.h"


vector<Node*> assignPoints(Gate* root, int n){
    vector<Node*> all_experts=root->getTerminalNodes();
    vector<Node*> z_final(n);
   
    for(int i=0;i<n;i++){
        int sub=rand() % all_experts.size();
        z_final[i]=all_experts[sub];
        cout<<"Point "<<i<<" assigned to "<<z_final[i]->name<<endl;
    }
    return z_final;
}

vec stdToArmaVec(vector<int> a){
    return conv_to<vec>::from(a);
}

void setUpData(mat A, vec* y, mat* X){
    *y=A.col(0);
    *X=A.submat(span(0,A.n_rows-1),span(1,A.n_cols-1));
    cout<<"Data loaded"<<endl;
}

double dnorm(double y, double mu, double sigma_sq){
     return exp(-0.5*(log(2*M_PI)+log(sigma_sq))-pow(y-mu,2)/(2*sigma_sq));
}

int main(){

cout<<"Set a random seed"<<endl;
arma_rng::set_seed(11);
//arma_rng::set_seed_random();

cout<<"Create a bunch of gates and experts"<<endl;
Gate* RP= new Gate(); 
Gate* G1=new Gate();
Gate* G2=new Gate();
Gate* G3=new Gate();
Gate* G4=new Gate();
Gate* G5=new Gate();
Expert* E1=new Expert();
Expert* E2=new Expert();
Expert* E3=new Expert();
Expert* E4=new Expert();
Expert* E5=new Expert();
Expert* E6=new Expert();
RP->name="RP";
G1->name="G1";
G2->name="G2";
G3->name="G3";
G4->name="G4";
G5->name="G5";
E1->name="E1";
E2->name="E2";
E3->name="E3";
E4->name="E4";
E5->name="E5";
E6->name="E6";

cout<<"Create a tree"<<endl;
RP->makeThisRootParent(G1);
G1->addChild(E1);
G1->addChild(G2);
G2->addChild(E2);
G2->addChild(E3);

cout<<"Issue IDs"<<endl;
G1->issueID();
G1->issueIDLR();

cout<<"Set expert models"<<endl;
NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;
E5->expertmodel=NF;
E6->expertmodel=NF;

cout<<"Load the data"<<endl;
mat A;
A.load("data_quad2.csv", csv_ascii);
vec y;
mat X;
setUpData(A,&y,&X);

cout<<"Set priors"<<endl;
vec mu_beta(X.n_cols); //prior mean for beta
mu_beta.zeros();
vec diags(X.n_cols);
diags.fill(0.5);
mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
mat Omega; //prior var-cov matrix for gamma, currently set in Gate.cpp, but should be taken outside
double a=1; //IG parameters for sigma^2
double b=2;
 
cout<<"Load initial allocations"<<endl;
vector<Node*> z_assign(y.size());
vec z_before(y.size());
z_before.load("z_before.csv", csv_ascii);

cout<<"Record point assignment as pointers"<<endl;
for(int i=0;i<z_before.size();i++){
     if(z_before[i]==1) z_assign[i]=E1;
     if(z_before[i]==2) z_assign[i]=E2;
     if(z_before[i]==3) z_assign[i]=E3;
}

cout<<"Estimate betas, set sigmas, estimate gammas"<<endl;
G1->estimateAllBetas(z_assign,y,X,0.1,mu_beta,Sigma_beta); 
cout<<"Betas estimated"<<endl;
G1->setAllSigmas(0.5);
cout<<"Sigmas set"<<endl;
G1->estimateAllGamas(z_assign,X,Omega);
cout<<"Gammas estimated"<<endl;


cout<<"Set up for MCMC"<<endl;
cout<<"Set the number of iterations"<<endl;
int N_MCMC=1000;
cout<<"Load new explanatory variable to predict on"<<endl;
mat X_new;
X_new.load("x_new.csv", csv_ascii);
int n_obs=X_new.n_rows;
cout<<"Set up a matrix to store predictions in"<<endl;
mat P(N_MCMC,n_obs);
cout<<"Set up matrices to store gammas in"<<endl;
mat gammas1(2,N_MCMC+1);
mat gammas2(2,N_MCMC+1);


cout<<"MCMC Start"<<endl;
vector<Node*> z_assign_new=G1->MCMC(N_MCMC,y,X,0,mu_beta,Sigma_beta,a,b,z_assign,&P,X_new,&gammas1,&gammas2);
cout<<"MCMC Finish"<<endl;

cout<<"Record parameters after MCMC"<<endl;
mat coefs(2,5);
coefs.col(0)=E1->beta;
coefs.col(1)=E2->beta;
coefs.col(2)=E3->beta;
coefs.col(3)=G1->gamma;
coefs.col(4)=G2->gamma;
mat sigmas(1,3);
sigmas.col(0)=E1->logsigma_sq;
sigmas.col(1)=E2->logsigma_sq;
sigmas.col(2)=E3->logsigma_sq;


cout<<"Record allocations after MCMC"<<endl;
vec z_afterMCMC(y.size());

for(int i=0;i<z_assign.size();i++){
    if(z_assign_new[i]->name=="E1") z_afterMCMC[i]=1;
    if(z_assign_new[i]->name=="E2") z_afterMCMC[i]=2;
    if(z_assign_new[i]->name=="E3") z_afterMCMC[i]=3;
}

cout<<"Save results"<<endl;
z_afterMCMC.save("z_afterMCMC.csv",csv_ascii);
P.save("P.csv",csv_ascii);
coefs.save("coefs.csv",csv_ascii);
sigmas.save("sigmas.csv",csv_ascii);
gammas1.save("gammas1.csv",csv_ascii);
gammas2.save("gammas2.csv",csv_ascii);
cout<<"Results saved"<<endl;

//G_1->swap(G_2,1,y2,X2);

vec x_record(2);
vec gamma_new(2);

z_assign_new=RP->split(y,X,E3,E4,G3,z_assign_new,0,0.5,0,50,&x_record,&gamma_new);
cout<<"Splitting done"<<endl;
cout<<"Record allocations after split"<<endl;
vec z_afterSplit(y.size());

for(int i=0;i<z_assign.size();i++){
    if(z_assign_new[i]->name=="E1") z_afterSplit[i]=1;
    if(z_assign_new[i]->name=="E2") z_afterSplit[i]=2;
    if(z_assign_new[i]->name=="E3") z_afterSplit[i]=3;
    if(z_assign_new[i]->name=="E4") z_afterSplit[i]=4;
}

mat coefs_new(2,7);
coefs_new.col(0)=E1->beta;
coefs_new.col(1)=E2->beta;
coefs_new.col(2)=E3->beta;
coefs_new.col(3)=E4->beta;
coefs_new.col(4)=G1->gamma;
coefs_new.col(5)=G2->gamma;
coefs_new.col(6)=G2->gamma;
mat sigmas_new(1,4);
sigmas_new.col(0)=E1->logsigma_sq;
sigmas_new.col(1)=E2->logsigma_sq;
sigmas_new.col(2)=E3->logsigma_sq;
sigmas_new.col(3)=E4->logsigma_sq;

z_afterSplit.save("z_afterSplit.csv",csv_ascii);
gamma_new.save("gamma_new.csv",csv_ascii);
x_record.save("x_record.csv",csv_ascii);
coefs_new.save("coefs_new.csv",csv_ascii);
sigmas_new.save("sigmas_new.csv",csv_ascii);
cout<<"Results saved"<<endl;

// dynamic_cast<Gate*>(RP->findNode("G1"))->printChildren();
// dynamic_cast<Gate*>(RP->findNode("G2"))->printChildren();
// dynamic_cast<Gate*>(RP->findNode("G3"))->printChildren();
// cout<<RP->findNode("E1")->Parent->name<<endl;
// cout<<RP->findNode("E2")->Parent->name<<endl;
// cout<<RP->findNode("E3")->Parent->name<<endl;
// cout<<RP->findNode("E4")->Parent->name<<endl;
// cout<<RP->findNode("G1")->Parent->name<<endl;
// cout<<RP->findNode("G2")->Parent->name<<endl;
// cout<<RP->findNode("G3")->Parent->name<<endl;

return 0;

}