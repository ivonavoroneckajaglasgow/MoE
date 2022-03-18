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
#include <chrono>

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


int main_mcycle(){ //mcycle data

Gate* RP=new Gate(); 
Gate* G1=new Gate();
Gate* G2=new Gate();
Gate* G3=new Gate();
Gate* G4=new Gate();
Expert* E1=new Expert();
Expert* E2=new Expert();
Expert* E3=new Expert();
Expert* E4=new Expert();
Expert* E5=new Expert();
RP->name="RP";
G1->name="G1";
G2->name="G2";
G3->name="G3";
G4->name="G4";
E1->name="E1";
E2->name="E2";
E3->name="E3";
E4->name="E4";
E5->name="E5";

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;
E5->expertmodel=NF;

mat A; vec y; mat X;
A.load("data_mcycle_stand.csv", csv_ascii);
setUpData(A,&y,&X);
vector<Node*> z_assign(y.size()); vec z_before(y.size());

//int Case=1;//start with 5 experts
int Case=2;//start with 1 expert

if(Case==1){
    RP->makeThisRootParent(G1);
    G1->addChild(G2);
    G1->addChild(G3);
    G2->addChild(E1);
    G2->addChild(E2);
    G3->addChild(E3);
    G3->addChild(G4);
    //G3->addChild(E4);
    G4->addChild(E4);
    G4->addChild(E5);
    // G1->addChild(E1);
    // G1->addChild(G2);
    // G2->addChild(E2);
    // G2->addChild(G3);
    // G3->addChild(E3);
    // G3->addChild(G4);
    // G4->addChild(E4);
    // G4->addChild(E5);
    z_before.load("z_mcycle.csv", csv_ascii);
}else{
    cout<<"Start with 1 expert"<<endl;
    RP->makeThisRootParent(G1);
    G1->addChild(E1);
    mat z_before(X.n_rows,1); z_before.fill(1);
    for(int i=0;i<(z_before.n_rows);i++) z_assign[i]=RP->findNode("E1");
}

G1->issueID();
G1->issueIDLR();

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(100); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
double gamma_var=50;
diags.fill(gamma_var); 
mat Omega=diagmat(diags); //prior var-cov matrix for gamma
Omega.submat(span(1,1),span(1,1))=100;
double a=1; double b=0.1; //prior IG for inverse gamma
//For the update of gamma, b refers to the rate

if(Case==1){
    for(int i=0;i<(z_before.n_rows);i++){
        if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
        if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
        if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
        if(as_scalar(z_before.row(i))==4) z_assign[i]=RP->findNode("E4");
        if(as_scalar(z_before.row(i))==5) z_assign[i]=RP->findNode("E5");
    }
}
     

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,-3,mu_beta,Sigma_beta); 
G1->setAllSigmas(0);
diags[0]=50,diags[1]=100;
mat Omega_helper=diagmat(diags);
if(Case==1){
    G1->estimateAllGamas(z_assign,X,Omega_helper);
}else{
    vec mu_gamma=mu_beta;
    G1->gamma= mu_gamma;
}

vector<Node*> gates=G1->getGates();
vector<Node*> experts=G1->getTerminalNodes();
for(int i=0;i<gates.size();i++) cout<<dynamic_cast<Gate*>(gates[i])->name<<dynamic_cast<Gate*>(gates[i])->gamma<<endl;
for(int i=0;i<experts.size();i++)cout<<dynamic_cast<Expert*>(experts[i])->name<<" beta: "<<dynamic_cast<Expert*>(experts[i])->beta<<"sigma: "<<dynamic_cast<Expert*>(experts[i])->logsigma_sq<<endl;

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(2*gamma_var); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=0.5;
dynamic_cast<Gate*>(RP->Children[0])->printChildren();

int N_MCMC=5000;
bool doRJ=1;
int RJ_every=10;
int L=1;
int predict_every=1;
int record_params_every=RJ_every;
mat accept_RJ(N_MCMC/RJ_every*L,2); accept_RJ.fill(5555);
//mat X_new=X;
mat X_new;
X_new.load("X_new_mcycle.csv", csv_ascii);
mat predictions(y.size(),N_MCMC/predict_every); predictions.fill(33333);
mat predictions_var(y.size(),N_MCMC/predict_every); predictions.fill(33333);
mat no_expt(N_MCMC/RJ_every*L,2); no_expt.fill(80000);
mat z_record(N_MCMC,y.size());

string penalty="poi"; //enter "poi" or "geo" is something else entered, no penalty applied
double lambda=pow(10,-1); //size prior

//Create matrices to store results in:
mat beta_record(X.n_cols+2,G1->countTerminals());beta_record.fill(-1);
for(int g=0;g<experts.size();g++){ beta_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Expert*>(experts[g])->beta;
}

mat sigma_record(3,G1->countTerminals()); sigma_record.fill(-1);
for(int g=0;g<experts.size();g++){ sigma_record.submat(span(1,1),span(g,g))=dynamic_cast<Expert*>(experts[g])->logsigma_sq;
}

mat gamma_record(X.n_cols+2,G1->countGates());gamma_record.fill(-1);
for(int g=0;g<gates.size();g++){gamma_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Gate*>(gates[g])->gamma;}

mat pi_record(2+X.n_rows,G1->countTerminals());pi_record.fill(-1);
for(int g=0;g<experts.size();g++){
    pi_record.submat(span(1,1),span(g,g))=experts[g]->id;
    pi_record.submat(span(2,pi_record.n_rows-1),span(g,g))=G1->getPathProb_mat(experts[g],X);
}
mat pi_record2;
pi_record2=pi_record;

mat dens(101,10);

cout<<"HERE IS THE BEGINNING"<<endl;
// Record start time
auto start = std::chrono::high_resolution_clock::now();
arma_rng::set_seed(3); //3, 321 and 87
vector<Node*> z_MCMC_RJ=G1->MCMC_RJ(N_MCMC, doRJ, RJ_every, L, &accept_RJ, y, X, mu_beta, Sigma_beta, a, b, Omega,mu_gamma1,Sigma_gamma1,sigma_epsilon,z_assign, X_new, &predictions,&predictions_var, predict_every, record_params_every, &no_expt,lambda,&z_record,&beta_record, &sigma_record, &gamma_record,&pi_record,&pi_record2,penalty,&dens);
// Record end time
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
cout<<"Complete"<<endl;
std::cout << "Elapsed time: " << elapsed.count() << " s\n";

accept_RJ.save("accept_RJ.csv",csv_ascii);
predictions.save("predictions.csv",csv_ascii);
predictions_var.save("predictions_var.csv",csv_ascii);
no_expt.save("no_expt.csv",csv_ascii);
z_record.save("z_record.csv",csv_ascii);
beta_record.save("beta_record.csv",csv_ascii);
sigma_record.save("sigma_record.csv",csv_ascii);
gamma_record.save("gamma_record.csv",csv_ascii);
pi_record.save("pi_record.csv",csv_ascii);
pi_record2.save("pi_record2.csv",csv_ascii);
dens.save("dens.csv",csv_ascii);
X_new.save("X_new.csv",csv_ascii);

return 0;

}

int main_sin(){ //going back to test a simple

Gate* RP=new Gate(); 
Gate* G1=new Gate();
Gate* G2=new Gate();
Gate* G3=new Gate();
Gate* G4=new Gate();
Expert* E1=new Expert();
Expert* E2=new Expert();
Expert* E3=new Expert();
Expert* E4=new Expert();
Expert* E5=new Expert();
RP->name="RP";
G1->name="G1";
G2->name="G2";
G3->name="G3";
G4->name="G4";
E1->name="E1";
E2->name="E2";
E3->name="E3";
E4->name="E4";
E5->name="E5";

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;
E5->expertmodel=NF;

mat A; vec y; mat X;
A.load("data_sin.csv", csv_ascii);
setUpData(A,&y,&X);
vector<Node*> z_assign(y.size()); vec z_before(y.size());

int Case=1;//srart with 3 experts

if(Case==1){
    RP->makeThisRootParent(G1);
    G1->addChild(E1);
    G1->addChild(G2);
    G2->addChild(E2);
    G2->addChild(E3);
    z_before.load("z_sin2.csv", csv_ascii);
}else{
    cout<<"else"<<endl;
}

G1->issueID();
G1->issueIDLR();

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(100); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
double gamma_var=50;
diags.fill(gamma_var); 
mat Omega=diagmat(diags); //prior var-cov matrix for gamma
double a=1; double b=10; //prior IG for inverse gamma

for(int i=0;i<(z_before.n_rows);i++){
      if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
      if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
      if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
      if(as_scalar(z_before.row(i))==4) z_assign[i]=RP->findNode("E4");
      if(as_scalar(z_before.row(i))==5) z_assign[i]=RP->findNode("E5");
}

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,-3,mu_beta,Sigma_beta); 
G1->setAllSigmas(-3);
G1->estimateAllGamas(z_assign,X,Omega);

vector<Node*> gates=G1->getGates();
vector<Node*> experts=G1->getTerminalNodes();
for(int i=0;i<gates.size();i++) cout<<dynamic_cast<Gate*>(gates[i])->name<<dynamic_cast<Gate*>(gates[i])->gamma<<endl;
for(int i=0;i<experts.size();i++) cout<<dynamic_cast<Expert*>(experts[i])->name<<dynamic_cast<Expert*>(experts[i])->beta<<endl;

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(2*gamma_var); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=0.5;
dynamic_cast<Gate*>(RP->Children[0])->printChildren();

int N_MCMC=5000;
bool doRJ=1;
int RJ_every=50;
int L=10;
int predict_every=1;
int record_params_every=50;
mat accept_RJ(N_MCMC/RJ_every*L,2); accept_RJ.fill(5555);
mat X_new=X; 
mat predictions(y.size(),N_MCMC/predict_every); predictions.fill(33333);
mat predictions_var(y.size(),N_MCMC/predict_every); predictions.fill(33333);
mat no_expt(N_MCMC/RJ_every*L,2); no_expt.fill(80000);
mat z_record(N_MCMC,y.size());

string penalty="poi"; //enter "poi" or "geo" is something else entered, no penalty applied
double lambda=pow(10,-1); //size prior

//Create matrices to store results in:
mat beta_record(X.n_cols+2,G1->countTerminals());beta_record.fill(-1);
for(int g=0;g<experts.size();g++){ beta_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Expert*>(experts[g])->beta;
}

mat sigma_record(3,G1->countTerminals()); sigma_record.fill(-1);
for(int g=0;g<experts.size();g++){ sigma_record.submat(span(1,1),span(g,g))=dynamic_cast<Expert*>(experts[g])->logsigma_sq;
}

mat gamma_record(X.n_cols+2,G1->countGates());gamma_record.fill(-1);
for(int g=0;g<gates.size();g++){gamma_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Gate*>(gates[g])->gamma;}

mat pi_record(2+X.n_rows,G1->countTerminals());pi_record.fill(-1);
for(int g=0;g<experts.size();g++){
    pi_record.submat(span(1,1),span(g,g))=experts[g]->id;
    pi_record.submat(span(2,pi_record.n_rows-1),span(g,g))=G1->getPathProb_mat(experts[g],X);
}
mat pi_record2;
pi_record2=pi_record;

mat dens(101,10);

// Record start time
auto start = std::chrono::high_resolution_clock::now();
arma_rng::set_seed(87); //3, 321 and 87
vector<Node*> z_MCMC_RJ=G1->MCMC_RJ(N_MCMC, doRJ, RJ_every, L, &accept_RJ, y, X, mu_beta, Sigma_beta, a, b, Omega,mu_gamma1,Sigma_gamma1,sigma_epsilon,z_assign, X_new, &predictions, &predictions_var, predict_every, record_params_every, &no_expt,lambda,&z_record,&beta_record, &sigma_record, &gamma_record,&pi_record,&pi_record2,penalty,&dens);
// Record end time
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
cout<<"Complete"<<endl;
std::cout << "Elapsed time: " << elapsed.count() << " s\n";

accept_RJ.save("accept_RJ.csv",csv_ascii);
predictions.save("predictions.csv",csv_ascii);
predictions_var.save("predictions_var.csv",csv_ascii);
no_expt.save("no_expt.csv",csv_ascii);
z_record.save("z_record.csv",csv_ascii);
beta_record.save("beta_record.csv",csv_ascii);
sigma_record.save("sigma_record.csv",csv_ascii);
gamma_record.save("gamma_record.csv",csv_ascii);
pi_record.save("pi_record.csv",csv_ascii);
pi_record2.save("pi_record2.csv",csv_ascii);
dens.save("dens.csv",csv_ascii);
//RJ_every.save("RJ_every.csv",csv_ascii);

// vec test(1000);a=0.1;b=10;
// for(int i=0;i<test.size();i++) test[i]=1/randg(distr_param(a,b));
// test.save("gamma_test.csv",csv_ascii);


return 0;

}

int main_real(){ //real estate

Gate* RP=new Gate(); 
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

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;
E5->expertmodel=NF;
E6->expertmodel=NF;

mat A; vec y; mat X;
A.load("RealEstateMRT.csv", csv_ascii);
setUpData(A,&y,&X);
vector<Node*> z_assign(y.size()); vec z_before(y.size());

int Case=1;//srart with 6 experts

if(Case==1){
    RP->makeThisRootParent(G1);
    G1->addChild(G2);
    G1->addChild(G3);
    G2->addChild(E1);
    G2->addChild(E2);
    G3->addChild(G4);
    G3->addChild(G5);
    G4->addChild(E3);
    G4->addChild(E4);
    G5->addChild(E5);
    G5->addChild(E6);
    z_before.load("z_RealEstateMRT.csv", csv_ascii);
}else{
    cout<<"else"<<endl;
}

G1->issueID();
G1->issueIDLR();

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(100); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
double gamma_var=50;
diags.fill(gamma_var); 
mat Omega=diagmat(diags); //prior var-cov matrix for gamma
Omega.submat(span(1,1),span(1,1))=100;
double a=1; double b=0.01; //prior IG for inverse gamma
//For the update of gamma, b refers to the rate

for(int i=0;i<(z_before.n_rows);i++){
      if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
      if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
      if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
      if(as_scalar(z_before.row(i))==4) z_assign[i]=RP->findNode("E4");
      if(as_scalar(z_before.row(i))==5) z_assign[i]=RP->findNode("E5");
      if(as_scalar(z_before.row(i))==6) z_assign[i]=RP->findNode("E6");
}

//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,1,mu_beta,Sigma_beta); 
G1->setAllSigmas(0);
diags[0]=50,diags[1]=100;
mat Omega_helper=diagmat(diags);
G1->estimateAllGamas(z_assign,X,Omega_helper);


vector<Node*> gates=G1->getGates();
vector<Node*> experts=G1->getTerminalNodes();
for(int i=0;i<gates.size();i++) cout<<dynamic_cast<Gate*>(gates[i])->name<<dynamic_cast<Gate*>(gates[i])->gamma<<endl;
for(int i=0;i<experts.size();i++)cout<<dynamic_cast<Expert*>(experts[i])->name<<" beta: "<<dynamic_cast<Expert*>(experts[i])->beta<<endl;//<<"sigma: "<<dynamic_cast<Expert*>(experts[i])->logsigma_sq<<endl;

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(2*gamma_var); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=0.5;
dynamic_cast<Gate*>(RP->Children[0])->printChildren();

int N_MCMC=5000;
bool doRJ=1;
int RJ_every=10;
int L=1;
int predict_every=1;
int record_params_every=RJ_every;
mat accept_RJ(N_MCMC/RJ_every*L,2); accept_RJ.fill(5555);
//mat X_new=X; 
mat X_new;
X_new.load("X_new_realestate.csv", csv_ascii);
mat predictions(y.size(),N_MCMC/predict_every); predictions.fill(33333);
mat predictions_var(y.size(),N_MCMC/predict_every); predictions.fill(33333);
mat no_expt(N_MCMC/RJ_every*L,2); no_expt.fill(80000);
mat z_record(N_MCMC,y.size());

string penalty="poi"; //enter "poi" or "geo" is something else entered, no penalty applied
double lambda=pow(10,-1); //size prior

//Create matrices to store results in:
mat beta_record(X.n_cols+2,G1->countTerminals());beta_record.fill(-1);
for(int g=0;g<experts.size();g++){ 
    beta_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Expert*>(experts[g])->beta;
    beta_record(3,g)=experts[g]->id;
}

mat sigma_record(3,G1->countTerminals()); sigma_record.fill(-1);
for(int g=0;g<experts.size();g++){ sigma_record.submat(span(1,1),span(g,g))=dynamic_cast<Expert*>(experts[g])->logsigma_sq;
}

mat gamma_record(X.n_cols+2,G1->countGates());gamma_record.fill(-1);
for(int g=0;g<gates.size();g++){gamma_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Gate*>(gates[g])->gamma;}

mat pi_record(2+X.n_rows,G1->countTerminals());pi_record.fill(-1);
for(int g=0;g<experts.size();g++){
    pi_record.submat(span(1,1),span(g,g))=experts[g]->id;
    pi_record.submat(span(2,pi_record.n_rows-1),span(g,g))=G1->getPathProb_mat(experts[g],X);
}
mat pi_record2;
pi_record2=pi_record;

mat dens(101,10);

cout<<"HERE IS THE BEGINNING"<<endl;
// Record start time
auto start = std::chrono::high_resolution_clock::now();
arma_rng::set_seed(3); //3, 321 and 87
vector<Node*> z_MCMC_RJ=G1->MCMC_RJ(N_MCMC, doRJ, RJ_every, L, &accept_RJ, y, X, mu_beta, Sigma_beta, a, b, Omega,mu_gamma1,Sigma_gamma1,sigma_epsilon,z_assign, X_new, &predictions, &predictions_var, predict_every, record_params_every, &no_expt,lambda,&z_record,&beta_record, &sigma_record, &gamma_record,&pi_record,&pi_record2,penalty,&dens);
// Record end time
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
cout<<"Complete"<<endl;
std::cout << "Elapsed time: " << elapsed.count() << " s\n";

dynamic_cast<Gate*>(RP->Children[0])->predictions_var(X_new);

accept_RJ.save("accept_RJ.csv",csv_ascii);
predictions.save("predictions.csv",csv_ascii);
predictions_var.save("predictions_var.csv",csv_ascii);
no_expt.save("no_expt.csv",csv_ascii);
z_record.save("z_record.csv",csv_ascii);
beta_record.save("beta_record.csv",csv_ascii);
sigma_record.save("sigma_record.csv",csv_ascii);
gamma_record.save("gamma_record.csv",csv_ascii);
pi_record.save("pi_record.csv",csv_ascii);
pi_record2.save("pi_record2.csv",csv_ascii);
dens.save("dens.csv",csv_ascii);
X_new.save("X_new.csv",csv_ascii);

// // vec test(100);
// // double c=1; double d=10;
// // arma_rng::set_seed(333); 
// // for(int i=0;i<100;i++) test[i]=1/randg( distr_param(c,d));
// // test.save("test.csv",csv_ascii);

// cout<<E1->expertmodel->IG_log(0.5,1,0.1)<<endl;

return 0;

}

int main(){ //real estate 2 vars

Gate* RP=new Gate(); 
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

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;
E4->expertmodel=NF;
E5->expertmodel=NF;
E6->expertmodel=NF;

mat A; vec y; mat X;
A.load("property2vars.csv", csv_ascii);
setUpData(A,&y,&X);
vector<Node*> z_assign(y.size()); vec z_before(y.size());

//int Case=1;//start with 6 experts
int Case=2;

if(Case==1){
     RP->makeThisRootParent(G1);
     G1->addChild(G2);
     G1->addChild(G3);
     G2->addChild(E1);
     G2->addChild(E2);
     G3->addChild(G4);
     G3->addChild(G5);
     G4->addChild(E3);
     G4->addChild(E4);
     G5->addChild(E5);
     G5->addChild(E6);
     z_before.load("z_RealEstateMRT.csv", csv_ascii);
 }else{
    cout<<"Start with 1 expert"<<endl;
    RP->makeThisRootParent(G1);
    G1->addChild(E1);
    mat z_before(X.n_rows,1); z_before.fill(1);
    for(int i=0;i<(z_before.n_rows);i++) z_assign[i]=RP->findNode("E1");
 }

 G1->issueID();
 G1->issueIDLR();

vec mu_beta(X.n_cols); mu_beta.zeros();
vec diags(X.n_cols); diags.fill(100); mat Sigma_beta=diagmat(diags); //prior var-cov matrix for beta
//mu_beta.print("mu_beta:"); Sigma_beta.print("Sigma_beta");
double gamma_var=50;
diags.fill(gamma_var); 
mat Omega=diagmat(diags); //prior var-cov matrix for gamma
//cout<<"Number of colums: "<<X.n_cols<<endl;
for(int i=0; i<(X.n_cols-1); i++) Omega.submat(span(i+1,i+1),span(i+1,i+1))=100;
//Omega.print("Omega:");
double a=1; double b=0.01; //prior IG for inverse gamma
//For the update of gamma, b refers to the rate

if(Case==1){
    for(int i=0;i<(z_before.n_rows);i++){
        if(as_scalar(z_before.row(i))==1) z_assign[i]=RP->findNode("E1");
        if(as_scalar(z_before.row(i))==2) z_assign[i]=RP->findNode("E2");
        if(as_scalar(z_before.row(i))==3) z_assign[i]=RP->findNode("E3");
        if(as_scalar(z_before.row(i))==4) z_assign[i]=RP->findNode("E4");
        if(as_scalar(z_before.row(i))==5) z_assign[i]=RP->findNode("E5");
        if(as_scalar(z_before.row(i))==6) z_assign[i]=RP->findNode("E6");
    }
}
//Set some initial values for parameters
G1->estimateAllBetas(z_assign,y,X,1,mu_beta,Sigma_beta); 
G1->setAllSigmas(0);
// diags[0]=50,diags[1]=100;
// mat Omega_helper=diagmat(diags);
if(Case==1){
    G1->estimateAllGamas(z_assign,X,Omega);
}else{
    vec mu_gamma=mu_beta;
    G1->gamma= mu_gamma;
}


vector<Node*> gates=G1->getGates();
vector<Node*> experts=G1->getTerminalNodes();
//for(int i=0;i<gates.size();i++) cout<<dynamic_cast<Gate*>(gates[i])->name<<dynamic_cast<Gate*>(gates[i])->gamma<<endl;
//for(int i=0;i<experts.size();i++)cout<<dynamic_cast<Expert*>(experts[i])->name<<" beta: "<<dynamic_cast<Expert*>(experts[i])->beta<<endl;//<<"sigma: "<<dynamic_cast<Expert*>(experts[i])->logsigma_sq<<endl;

//Parameters for forward jump proposal
vec mu_gamma1(X.n_cols-1);mu_gamma1.fill(0);
vec diags2(X.n_cols-1); diags2.fill(2*gamma_var); mat Sigma_gamma1=diagmat(diags2); 
double sigma_epsilon=0.5;
//mu_gamma1.print("mu_gamma1: ");Sigma_gamma1.print("Sigma_gamma1: ");

int N_MCMC=2000;
bool doRJ=1;
int RJ_every=10;
int L=1;
int predict_every=1;
int record_params_every=RJ_every;
mat accept_RJ(N_MCMC/RJ_every*L,2); accept_RJ.fill(5555);
//mat X_new=X; 
mat X_new;
X_new.load("X_new_realestate_2vars.csv", csv_ascii);
mat predictions(X_new.n_rows,N_MCMC/predict_every); predictions.fill(33333);
mat predictions_var(X_new.n_rows,N_MCMC/predict_every); predictions.fill(33333);
mat no_expt(N_MCMC/RJ_every*L,2); no_expt.fill(80000);
mat z_record(N_MCMC,y.size());

string penalty="poi"; //enter "poi" or "geo" is something else entered, no penalty applied
double lambda=pow(10,-1); //size prior

//Create matrices to store results in:
mat beta_record(X.n_cols+2,G1->countTerminals());beta_record.fill(-1);
for(int g=0;g<experts.size();g++){ 
     beta_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Expert*>(experts[g])->beta;
     beta_record(beta_record.n_rows-1,g)=experts[g]->id;
}
//beta_record.print("Starter beta_record:");
mat sigma_record(3,G1->countTerminals()); sigma_record.fill(-1);
for(int g=0;g<experts.size();g++){ sigma_record.submat(span(1,1),span(g,g))=dynamic_cast<Expert*>(experts[g])->logsigma_sq;
}
//sigma_record.print("Starter sigma_record:");

mat gamma_record(X.n_cols+2,G1->countGates());gamma_record.fill(-1);
for(int g=0;g<gates.size();g++){gamma_record.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Gate*>(gates[g])->gamma;}
//gamma_record.print("Starter gamma_record:");

mat pi_record(2+X.n_rows,G1->countTerminals());pi_record.fill(-1);
for(int g=0;g<experts.size();g++){
    pi_record.submat(span(1,1),span(g,g))=experts[g]->id;
    pi_record.submat(span(2,pi_record.n_rows-1),span(g,g))=G1->getPathProb_mat(experts[g],X);
}
mat pi_record2;
pi_record2=pi_record;
//pi_record.print("Starter pi_record:");

mat dens(101,10);

cout<<"HERE IS THE BEGINNING"<<endl;
// Record start time
auto start = std::chrono::high_resolution_clock::now();
arma_rng::set_seed(3); //3, 321 and 87
vector<Node*> z_MCMC_RJ=G1->MCMC_RJ(N_MCMC, doRJ, RJ_every, L, &accept_RJ, y, X, mu_beta, Sigma_beta, a, b, Omega,mu_gamma1,Sigma_gamma1,sigma_epsilon,z_assign, X_new, &predictions, &predictions_var, predict_every, record_params_every, &no_expt,lambda,&z_record,&beta_record, &sigma_record, &gamma_record,&pi_record,&pi_record2,penalty,&dens);
// Record end time
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
cout<<"Complete"<<endl;
std::cout << "Elapsed time: " << elapsed.count() << " s\n";

accept_RJ.save("accept_RJ.csv",csv_ascii);
predictions.save("predictions.csv",csv_ascii);
predictions_var.save("predictions_var.csv",csv_ascii);
no_expt.save("no_expt.csv",csv_ascii);
z_record.save("z_record.csv",csv_ascii);
beta_record.save("beta_record.csv",csv_ascii);
sigma_record.save("sigma_record.csv",csv_ascii);
gamma_record.save("gamma_record.csv",csv_ascii);
pi_record.save("pi_record.csv",csv_ascii);
pi_record2.save("pi_record2.csv",csv_ascii);
dens.save("dens.csv",csv_ascii);
X_new.save("X_new.csv",csv_ascii);


return 0;

}

