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
z_afterMCMC.save("z_afterMCMC_split1.csv",csv_ascii);

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

int main_split_testing(){
    int N_split=1000; vec split_record(N_split); split_record.fill(0);
    mat beta1(2,1); mat beta2(2,1); mat beta3(2,1); 
    beta1.fill(-100);beta2.fill(-100);beta3.fill(-100);
    mat gamma1(2,1); mat gamma2(2,1); mat z(100,1);
    gamma1.fill(-100); gamma2.fill(-100); z.fill(-100);

    mat beta1MCMC(2,1); mat beta2MCMC(2,1); mat beta3MCMC(2,1); 
    beta1MCMC.fill(-100);beta2MCMC.fill(-100);beta3MCMC.fill(-100);
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

int one_merge_proposal(){
    
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
G1->addChild(G2);
G2->addChild(E2);
G2->addChild(E3);

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
z_before.load("z_for_merge.csv", csv_ascii);
for(int i=0;i<(z_before.n_rows);i++){
    if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
    if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
    if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
}

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,0.1,mu_beta,Sigma_beta); 
G1->setAllSigmas(0.5);
G1->estimateAllGamas(z_assign,X,Omega);

//Set up for MCMC
int N_MCMC=10;
bool doMCMC=1;
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
      if(z_assign_new[i]->name=="E3") z_afterMCMC[i]=3;
}
z_afterMCMC.save("z_afterMCMC_merge1.csv",csv_ascii);

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(20); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=1;

vec z_afterMerge(y.size());
//arma_rng::set_seed_random();
vector<Node*> z_post_merge=G1->merge(y,X,G2,z_assign_new,mu_beta,Sigma_beta,a,b,Omega);
    
for(int i=0;i<z_assign.size();i++){
    if(z_post_merge[i]->name=="E1") z_afterMerge[i]=1;
    if(z_post_merge[i]->name=="E2") z_afterMerge[i]=2;
    if(z_post_merge[i]->name=="E3") z_afterMerge[i]=3;
}
//z_afterSplit.save("z_afterSplit.csv",csv_ascii);

    if(any(z_afterMerge==3)){
        cout<<"Merge rejected."<<endl;
        return 0;
    }else{
        cout<<"Merge accepted."<<endl;
        return 1;
    }
    //return 0;
}

int main_merge_testing(){ //Merge, where it should be rejected
    int N_merge=1000; vec merge_record(N_merge); merge_record.fill(0);

    for(int k=0;k<N_merge;k++){
        cout<<"Proposing merge number " <<k<<endl;
        merge_record[k]=one_merge_proposal();
    }
    merge_record.save("merge_record.csv",csv_ascii);

return 0;
}

int one_merge_proposal2(){
    
Gate* RP=new Gate(); 
Gate* G1=new Gate();
Gate* G2=new Gate();
Gate* G3=new Gate();
Expert* E1=new Expert();
Expert* E2=new Expert();
Expert* E3=new Expert();
Expert* E4=new Expert();
RP->name="RP";
G1->name="G1";
G2->name="G2";
G3->name="G3";
E1->name="E1";
E2->name="E2";
E3->name="E3";
E4->name="E4";

RP->makeThisRootParent(G1);
G1->addChild(G2);
G1->addChild(G3);
G2->addChild(E1);
G2->addChild(E2);
G3->addChild(E3);
G3->addChild(E4);

G1->issueID();
G1->issueIDLR();

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;

mat A; vec y; mat X;
A.load("data_sin_stand.csv", csv_ascii);
setUpData(A,&y,&X);

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(0.5); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
diags.fill(1); mat Omega=diagmat(diags); //prior var-cov matrix for gamma
double a=0.01; double b=0.01; //prior IG for inverse gamma


vector<Node*> z_assign(y.size()); vec z_before(y.size());
z_before.load("z_4exp_merge.csv", csv_ascii);
for(int i=0;i<(z_before.n_rows);i++){
    if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
    if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
    if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
    if(as_scalar(z_before.row(i))==4) z_assign[i]=RP->findNode("E4");
}

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,0.1,mu_beta,Sigma_beta); 
G1->setAllSigmas(0.5);
G1->estimateAllGamas(z_assign,X,Omega);

//Set up for MCMC
int N_MCMC=10;
bool doMCMC=1;
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
      if(z_assign_new[i]->name=="E3") z_afterMCMC[i]=3;
      if(z_assign_new[i]->name=="E4") z_afterMCMC[i]=4;
}
z_afterMCMC.save("z_afterMCMC_merge2.csv",csv_ascii);

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(20); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=1;

vec z_afterMerge(y.size());
//arma_rng::set_seed_random();
vector<Node*> z_post_merge=G1->merge(y,X,G3,z_assign_new,mu_beta,Sigma_beta,a,b,Omega);
    
for(int i=0;i<z_assign.size();i++){
    if(z_post_merge[i]->name=="E1") z_afterMerge[i]=1;
    if(z_post_merge[i]->name=="E2") z_afterMerge[i]=2;
    if(z_post_merge[i]->name=="E3") z_afterMerge[i]=3;
    if(z_post_merge[i]->name=="E4") z_afterMerge[i]=4;
}
//z_afterSplit.save("z_afterSplit.csv",csv_ascii);

    if(any(z_afterMerge==4)){
        cout<<"Merge rejected."<<endl;
        return 0;
    }else{
        cout<<"Merge accepted."<<endl;
        return 1;
    }
    //return 0;
}

int main_merge_testing2(){ //merge where it should be accepted
    int N_merge=1000; vec merge_record(N_merge); merge_record.fill(0);

    for(int k=0;k<N_merge;k++){
        cout<<"Proposing merge number " <<k<<endl;
        merge_record[k]=one_merge_proposal2();
    }
    merge_record.save("merge_record2.csv",csv_ascii);

return 0;
}

int one_split_proposal2(mat* beta1, mat* beta2, mat* beta3, mat* beta4, mat* gamma1, mat* gamma2, mat* gamma3, mat* z, 
                        mat* beta1MCMC, mat* beta2MCMC, mat* beta3MCMC, mat* beta4MCMC, mat* gamma1MCMC, mat* gamma2MCMC, mat* gamma3MCMC, mat* zMCMC){
    
Gate* RP=new Gate(); 
Gate* G1=new Gate();
Gate* G2=new Gate();
Gate* G3=new Gate();
Expert* E1=new Expert();
Expert* E2=new Expert();
Expert* E3=new Expert();
Expert* E4=new Expert();
RP->name="RP";
G1->name="G1";
G2->name="G2";
G3->name="G3";
E1->name="E1";
E2->name="E2";
E3->name="E3";
E4->name="E4";

RP->makeThisRootParent(G1);
G1->addChild(E1);
G1->addChild(G2);
G2->addChild(E2);
G2->addChild(E3);

G1->issueID();
G1->issueIDLR();

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;

mat A; vec y; mat X;
A.load("data_sin_stand.csv", csv_ascii);
setUpData(A,&y,&X);

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(0.5); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
diags.fill(1); mat Omega=diagmat(diags); //prior var-cov matrix for gamma
double a=0.01; double b=0.01; //prior IG for inverse gamma


vector<Node*> z_assign(y.size()); vec z_before(y.size());
z_before.load("z_expected.csv", csv_ascii);
for(int i=0;i<(z_before.n_rows);i++){
      if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
      if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
      if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
}

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,0.1,mu_beta,Sigma_beta); 
G1->setAllSigmas(0.5);
G1->estimateAllGamas(z_assign,X,Omega);

//Set up for MCMC
int N_MCMC=10;
bool doMCMC=1;
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
      if(z_assign_new[i]->name=="E3") z_afterMCMC[i]=3;
}
z_afterMCMC.save("z_afterMCMC_split2.csv",csv_ascii);

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(20); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=1;

vec z_afterSplit(y.size());
//arma_rng::set_seed_random();
vector<Node*> z_post_split=G1->split2(y,X,E3,E4,G3,z_assign_new,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
    
for(int i=0;i<z_assign.size();i++){
    if(z_post_split[i]->name=="E1") z_afterSplit[i]=1;
    if(z_post_split[i]->name=="E2") z_afterSplit[i]=2;
    if(z_post_split[i]->name=="E3") z_afterSplit[i]=3;
    if(z_post_split[i]->name=="E4") z_afterSplit[i]=4;
}
//z_afterSplit.save("z_afterSplit.csv",csv_ascii);

    if(any(z_afterSplit==4)){
   
        (*beta1).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E1"))->beta);
        (*beta2).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E2"))->beta);
        (*beta3).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E3"))->beta);
        (*beta4).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E4"))->beta);
        (*gamma1).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G1"))->gamma);
        (*gamma2).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G2"))->gamma);
        (*gamma3).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G3"))->gamma);
        (*z).insert_cols(0,z_afterSplit);
 
        cout<<"Performing MCMC after accepted jump."<<endl;
        vector<Node*> z_post_MCMC=G1->MCMC(100,y,X,mu_beta,Sigma_beta,a,b,Omega,z_assign);
        cout<<"Done."<<endl;
        vec z_afterMCMC(y.size());
       
        for(int i=0;i<z_assign.size();i++){
            if(z_post_MCMC[i]->name=="E1") z_afterMCMC[i]=1;
            if(z_post_MCMC[i]->name=="E2") z_afterMCMC[i]=2;
            if(z_post_MCMC[i]->name=="E3") z_afterMCMC[i]=3;
            if(z_post_MCMC[i]->name=="E4") z_afterMCMC[i]=4;
        }

        (*beta1MCMC).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E1"))->beta);
        (*beta2MCMC).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E2"))->beta);
        (*beta3MCMC).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E3"))->beta);
        (*beta4MCMC).insert_cols(0,dynamic_cast<Expert*>(RP->findNode("E3"))->beta);
        (*gamma1MCMC).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G1"))->gamma);
        (*gamma2MCMC).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G2"))->gamma);
        (*gamma3MCMC).insert_cols(0,dynamic_cast<Gate*>(RP->findNode("G3"))->gamma);
        (*zMCMC).insert_cols(0,z_afterMCMC);

        return 1;
    }else{
        //cout<<"Split rejected"<<endl;
        return 0;
    }
}

int main_split_test(){
    int N_split=1000; vec split_record(N_split); split_record.fill(0);
    mat beta1(2,1); mat beta2(2,1); mat beta3(2,1); mat beta4(2,1); 
    beta1.fill(-100);beta2.fill(-100);beta3.fill(-100);beta4.fill(-100);
    mat gamma1(2,1); mat gamma2(2,1); mat gamma3(2,1); mat z(100,1);
    gamma1.fill(-100); gamma2.fill(-100); gamma3.fill(-100); z.fill(-100);


    mat beta1MCMC(2,1); mat beta2MCMC(2,1); mat beta3MCMC(2,1); mat beta4MCMC(2,1);
    beta1MCMC.fill(-100);beta2MCMC.fill(-100);beta3MCMC.fill(-100);beta4MCMC.fill(-100);
    mat gamma1MCMC(2,1); mat gamma2MCMC(2,1); mat gamma3MCMC(2,1); mat zMCMC(100,1);
    gamma1MCMC.fill(-100); gamma2MCMC.fill(-100); gamma3MCMC.fill(-100); zMCMC.fill(-100);

    for(int k=0;k<N_split;k++){
        cout<<"Proposing forward jump number " <<k<<endl;
        split_record[k]=one_split_proposal2(&beta1,&beta2,&beta3,&beta4,&gamma1,&gamma2,&gamma3,&z,&beta1MCMC,&beta2MCMC,&beta3MCMC,&beta4MCMC,&gamma1MCMC,&gamma2MCMC,&gamma3MCMC,&zMCMC);
    }
    split_record.save("split_record2.csv",csv_ascii);
    beta1.save("beta12.csv",csv_ascii);
    beta2.save("beta22.csv",csv_ascii);
    beta3.save("beta32.csv",csv_ascii);
    beta4.save("beta42.csv",csv_ascii);
    gamma1.save("gamma12.csv",csv_ascii);
    gamma2.save("gamma22.csv",csv_ascii);
    gamma3.save("gamma32.csv",csv_ascii);
    z.save("z2.csv",csv_ascii);
    beta1MCMC.save("beta1MCMC2.csv",csv_ascii);
    beta2MCMC.save("beta2MCMC2.csv",csv_ascii);
    beta3MCMC.save("beta3MCMC2.csv",csv_ascii);
    beta4MCMC.save("beta4MCMC2.csv",csv_ascii);
    gamma1MCMC.save("gamma1MCMC2.csv",csv_ascii);
    gamma2MCMC.save("gamma2MCMC2.csv",csv_ascii);
    gamma3MCMC.save("gamma3MCMC2.csv",csv_ascii);
    zMCMC.save("zMCMC2.csv",csv_ascii);
}

int main_test(){
mat A; vec y; mat X;
A.load("data_sin_stand.csv", csv_ascii);
setUpData(A,&y,&X);

Expert* E1=new Expert();
Expert* E2=new Expert();
NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
Gate* RP=new Gate(); 
Gate* G1=new Gate();
RP->makeThisRootParent(G1);
G1->addChild(E1);
G1->addChild(E2);
G1->issueID();
G1->issueIDLR();

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(0.5); mat Sigma_beta=diagmat(diags);
diags.fill(1); mat Omega=diagmat(diags);
vec beta("0.5, 1");
double logsigma_sq=log(0.5);
double a=0.01;
double b=0.01;

G1->gamma=beta;

vector<Node*> z_assign(y.size());
for(int i=0; i<z_assign.size();i++) z_assign[i]=E1;

//cout<<E1->expertmodel->qSigma(y, X, beta, logsigma_sq, a, b)<<endl;
//cout<<E1->expertmodel->qBeta(y, X, beta, logsigma_sq, mu_beta, Sigma_beta)<<endl;
//cout<<G1->qGamma(X,z_assign, Omega)<<endl;

return 0;
}

int main(){ //testing forward and backward together
int Case=3;

Gate* RP=new Gate(); 
Gate* G1=new Gate();
Gate* G2=new Gate();
Gate* G3=new Gate();
Expert* E1=new Expert();
Expert* E2=new Expert();
Expert* E3=new Expert();
Expert* E4=new Expert();
RP->name="RP";
G1->name="G1";
G2->name="G2";
G3->name="G3";
E1->name="E1";
E2->name="E2";
E3->name="E3";
E4->name="E4";

mat A; vec y; mat X;
A.load("data_sin_stand.csv", csv_ascii);
setUpData(A,&y,&X);
vector<Node*> z_assign(y.size()); vec z_before(y.size());

if(Case==1){
//Case 1) Start 2, Desired 3
    RP->makeThisRootParent(G1);
    G1->addChild(E1);
    G1->addChild(E2);
    z_before.load("z_before.csv", csv_ascii); //Case 1
}else{
    if(Case==2){
//Case 2) Start 3, Desired 3
        RP->makeThisRootParent(G1);
        G1->addChild(E1);
        G1->addChild(G2);
        G2->addChild(E2);
        G2->addChild(E3);
        z_before.load("z_for_merge.csv", csv_ascii); //Case 2
    }else{
//Case 3) Start 4, Desired 3
        RP->makeThisRootParent(G1);
        G1->addChild(G2);
        G1->addChild(G3);
        G2->addChild(E1);
        G2->addChild(E2);
        G3->addChild(E3);
        G3->addChild(E4);
        z_before.load("z_4exp_merge.csv", csv_ascii); //Case 3      
    }
}

G1->issueID();
G1->issueIDLR();

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;


vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(0.5); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
diags.fill(1); mat Omega=diagmat(diags); //prior var-cov matrix for gamma
double a=0.01; double b=0.01; //prior IG for inverse gamma


for(int i=0;i<(z_before.n_rows);i++){
      if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
      if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
      if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
      if(as_scalar(z_before.row(i))==4) z_assign[i]=RP->findNode("E4");
}

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,0.1,mu_beta,Sigma_beta); 
G1->setAllSigmas(0.5);
G1->estimateAllGamas(z_assign,X,Omega);

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(20); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=1;

arma_rng::set_seed(3);
vector<Node*> z_MCMC_RJ=G1->MCMC_RJ(5000,500, y, X, mu_beta, Sigma_beta, a, b, Omega,mu_gamma1,Sigma_gamma1,sigma_epsilon,z_assign);
cout<<"Complete"<<endl;

vec z_MCMCRJ(y.size());
for(int i=0;i<y.size();i++){
    if(z_MCMC_RJ[i]->name=="E1") z_MCMCRJ[i]=1;
    if(z_MCMC_RJ[i]->name=="E2") z_MCMCRJ[i]=2;
    if(z_MCMC_RJ[i]->name=="E3") z_MCMCRJ[i]=3;
    if(z_MCMC_RJ[i]->name=="E4") z_MCMCRJ[i]=4;
    if(z_MCMC_RJ[i]->name=="E5") z_MCMCRJ[i]=5;
    if(z_MCMC_RJ[i]->name=="E6") z_MCMCRJ[i]=6;
    if(z_MCMC_RJ[i]->name=="E7") z_MCMCRJ[i]=7;
    if(z_MCMC_RJ[i]->name=="E8") z_MCMCRJ[i]=8;
    if(z_MCMC_RJ[i]->name=="E9") z_MCMCRJ[i]=9;
    if(z_MCMC_RJ[i]->name=="E10") z_MCMCRJ[i]=10;
}

z_MCMCRJ.save("z_MCMCRJ.csv",csv_ascii);

    return 0;
}