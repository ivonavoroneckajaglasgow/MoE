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
    //cout<<"Data loaded"<<endl;
}

double dnorm(double y, double mu, double sigma_sq){
     return exp(-0.5*(log(2*M_PI)+log(sigma_sq))-pow(y-mu,2)/(2*sigma_sq));
}


int one_split_proposal(mat* beta1, mat* beta2, mat* beta3, mat* gamma1, mat* gamma2, mat* z, 
                       mat* beta1MCMC, mat* beta2MCMC, mat* beta3MCMC, mat* gamma1MCMC, mat* gamma2MCMC, mat* zMCMC){
    
Gate* RP=new Gate(); 
Gate* G1=new Gate();
Gate* G2=new Gate();
Expert* E1=new Expert();
Expert* E2=new Expert();
Expert* E3=new Expert();
RP->name="RP";
G1->name="G1";
G2->name="G2";
E1->name="E1";
E2->name="E2";
E3->name="E3";

RP->makeThisRootParent(G1);
G1->addChild(E1);
G1->addChild(E2);

G1->issueID();
G1->issueIDLR();

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;

mat A; vec y; mat X;
A.load("data_sin_stand.csv", csv_ascii);
setUpData(A,&y,&X);

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(0.5); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
diags.fill(1); mat Omega=diagmat(diags); //prior var-cov matrix for gamma
double a=0.01; double b=0.01; //prior IG for inverse gamma


vector<Node*> z_assign(y.size()); vec z_before(y.size());
z_before.load("z_before.csv", csv_ascii);
for(int i=0;i<(z_before.n_rows);i++){
      if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
      if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
}

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,0.1,mu_beta,Sigma_beta); 
G1->setAllSigmas(0.5);
G1->estimateAllGamas(z_assign,X,Omega);

//Set up for MCMC
int N_MCMC=10;
bool doMCMC=0;
vector<Node*> z_assign_new; 
if(doMCMC==1){
    z_assign_new=G1->MCMC(N_MCMC,y,X,mu_beta,Sigma_beta,a,b,Omega,z_assign);
}else{
    z_assign_new=z_assign;
}

vec z_afterMCMC(y.size());
for(int i=0;i<z_assign.size();i++){
      if(z_assign_new[i]->name=="E1") z_afterMCMC[i]=1;
      if(z_assign_new[i]->name=="E2") z_afterMCMC[i]=2;
}
z_afterMCMC.save("z_afterMCMC.csv",csv_ascii);

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(20); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=1;

vec z_afterSplit(y.size());
//arma_rng::set_seed_random();
vector<Node*> z_post_split=G1->split2(y,X,E2,E3,G2,z_assign_new,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
    
for(int i=0;i<z_assign.size();i++){
    if(z_post_split[i]->name=="E1") z_afterSplit[i]=1;
    if(z_post_split[i]->name=="E2") z_afterSplit[i]=2;
    if(z_post_split[i]->name=="E3") z_afterSplit[i]=3;
}
//z_afterSplit.save("z_afterSplit.csv",csv_ascii);

    if(any(z_afterSplit==3)){
        (*beta1).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E1"))->beta);
        (*beta2).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E2"))->beta);
        (*beta3).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E3"))->beta);
        (*gamma1).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G1"))->gamma);
        (*gamma2).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G2"))->gamma);
        (*z).insert_cols(0,z_afterSplit);

        cout<<"Performing MCMC after accepted jump."<<endl;
        vector<Node*> z_post_MCMC=G1->MCMC(100,y,X,mu_beta,Sigma_beta,a,b,Omega,z_assign);
        cout<<"Done."<<endl;
        vec z_afterMCMC(y.size());
       
        for(int i=0;i<z_assign.size();i++){
            if(z_post_MCMC[i]->name=="E1") z_afterMCMC[i]=1;
            if(z_post_MCMC[i]->name=="E2") z_afterMCMC[i]=2;
            if(z_post_MCMC[i]->name=="E3") z_afterMCMC[i]=3;
        }

        (*beta1MCMC).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E1"))->beta);
        (*beta2MCMC).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E2"))->beta);
        (*beta3MCMC).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E3"))->beta);
        (*gamma1MCMC).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G1"))->gamma);
        (*gamma2MCMC).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G2"))->gamma);
        (*zMCMC).insert_cols(0,z_afterMCMC);

        return 1;
    }else{
        return 0;
    }
}

int main(){
    int N_split=1000; vec split_record(N_split); split_record.fill(0);
    mat beta1(2,1); mat beta2(2,1); mat beta3(2,1); 
    beta1.fill(-100);beta2.fill(-100);beta2.fill(-100);
    mat gamma1(2,1); mat gamma2(2,1); mat z(100,1);
    gamma1.fill(-100); gamma2.fill(-100); z.fill(-100);

    mat beta1MCMC(2,1); mat beta2MCMC(2,1); mat beta3MCMC(2,1); 
    beta1MCMC.fill(-100);beta2MCMC.fill(-100);beta2MCMC.fill(-100);
    mat gamma1MCMC(2,1); mat gamma2MCMC(2,1); mat zMCMC(100,1);
    gamma1MCMC.fill(-100); gamma2MCMC.fill(-100); zMCMC.fill(-100);

    for(int k=0;k<N_split;k++){
        cout<<"Proposing forward jump number " <<k<<endl;
        split_record[k]=one_split_proposal(&beta1,&beta2,&beta3,&gamma1,&gamma2,&z,&beta1MCMC,&beta2MCMC,&beta3MCMC,&gamma1MCMC,&gamma2MCMC,&zMCMC);
    }
    split_record.save("split_record.csv",csv_ascii);
    beta1.save("beta1.csv",csv_ascii);
    beta2.save("beta2.csv",csv_ascii);
    beta3.save("beta3.csv",csv_ascii);
    gamma1.save("gamma1.csv",csv_ascii);
    gamma2.save("gamma2.csv",csv_ascii);
    z.save("z.csv",csv_ascii);
    beta1MCMC.save("beta1MCMC.csv",csv_ascii);
    beta2MCMC.save("beta2MCMC.csv",csv_ascii);
    beta3MCMC.save("beta3MCMC.csv",csv_ascii);
    gamma1MCMC.save("gamma1MCMC.csv",csv_ascii);
    gamma2MCMC.save("gamma2MCMC.csv",csv_ascii);
    zMCMC.save("zMCMC.csv",csv_ascii);
}