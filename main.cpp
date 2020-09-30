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

int main(){

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

//First Tree
G1->addChild(G2);
G1->addChild(E1);
G2->addChild(E2);
G2->addChild(E3);

G1->issueIDLR();

//Second Tree
G3->addChild(E4);
G3->addChild(G4);
G4->addChild(E5);
G4->addChild(E6);

G3->issueIDLR();

int n=10;
int p=2;

mat X(n,p,fill::randn);
vec y(n, fill::randn);

cout<<"Next, create some data with "<< n<<" points and "<<p<<" explanatory variables."<<endl;
cout<<"Now randomly assign the "<<n <<" points to the experts in the tree."<<endl;

vector<Node*> z_tree1=assignPoints(G1,n);
vector<Node*> z_tree2=assignPoints(G3,n);

cout<<"Now, I will estimate gammas for the gates based on the point allocation:"<<endl;

mat z_G1=G1->getZ(z_tree1);
mat z_G2=G2->getZ(z_tree1);
mat z_G3=G3->getZ(z_tree2);
mat z_G4=G4->getZ(z_tree2);

mat Omega;

cout<<"Now estimate gammas"<<endl;
G1->gamma=G1->findGammaMLE(G1->subsetX(X,G1->getPointIndices(z_tree1)),z_G1,Omega);
G2->gamma=G2->findGammaMLE(G2->subsetX(X,G2->getPointIndices(z_tree1)),z_G2,Omega);
G3->gamma=G3->findGammaMLE(G3->subsetX(X,G3->getPointIndices(z_tree2)),z_G3,Omega);
G4->gamma=G4->findGammaMLE(G4->subsetX(X,G4->getPointIndices(z_tree2)),z_G4,Omega);

//Set tree 1 to be Normal Family
//Set tree 2 to be Normal Model

NormalFamily* NF=new NormalFamily();
E1->expertmodel=NF;
E2->expertmodel=NF;
E3->expertmodel=NF;

NormalModel* NM=new NormalModel();
E4->expertmodel=NM;
E5->expertmodel=NM;
E6->expertmodel=NM;

double logsigma_sq=1;
vec mu_beta(X.n_cols);
mu_beta.zeros();
vec diags(X.n_cols);
diags.fill(0.0001);
mat Sigma_beta=diagmat(diags);

cout<<"Next, estimate betas for Normal Family experts based on the point allocations."<<endl;
E1->beta=E1->expertmodel->findBeta(E1->subsetY(y,E1->getPointIndices(z_tree1)),E1->subsetX(X,E1->getPointIndices(z_tree1)),logsigma_sq,mu_beta,Sigma_beta);
E2->beta=E2->expertmodel->findBeta(E2->subsetY(y,E2->getPointIndices(z_tree1)),E2->subsetX(X,E2->getPointIndices(z_tree1)),logsigma_sq,mu_beta,Sigma_beta);
E3->beta=E3->expertmodel->findBeta(E3->subsetY(y,E3->getPointIndices(z_tree1)),E3->subsetX(X,E3->getPointIndices(z_tree1)),logsigma_sq,mu_beta,Sigma_beta);

cout<<"Now estimate betas for Normal Model using ML."<<endl;
E4->beta=E4->expertmodel->findBeta(E4->subsetY(y,E4->getPointIndices(z_tree2)),E4->subsetX(X,E4->getPointIndices(z_tree2)),logsigma_sq);
E5->beta=E5->expertmodel->findBeta(E5->subsetY(y,E5->getPointIndices(z_tree2)),E5->subsetX(X,E5->getPointIndices(z_tree2)),logsigma_sq);
E6->beta=E6->expertmodel->findBeta(E6->subsetY(y,E6->getPointIndices(z_tree2)),E6->subsetX(X,E6->getPointIndices(z_tree2)),logsigma_sq);

cout<<"For simplicity, keep all logsigmasq="<<1<<" for all experts."<<endl;

E1->logsigma_sq=1;
E2->logsigma_sq=1;
E3->logsigma_sq=1;
E4->logsigma_sq=1;
E5->logsigma_sq=1;
E6->logsigma_sq=1;

double a=1;
double b=2;

int N=10;

cout<<"Try MCMC for tree 1."<<endl;
vector<Node*> z_tree1_updated=G1->MCMC(N,y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_tree1);
for(int i=0;i<z_tree1.size();i++){
    cout<<"Point "<<i<<" initially was in "<<z_tree1[i]->name<<" and now is in "<<z_tree1_updated[i]->name<<endl;
}

cout<<"Try MCMC for tree 2."<<endl;
vector<Node*> z_tree2_updated=G3->MCMC(N,y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_tree1);
for(int i=0;i<z_tree2.size();i++){
    cout<<"Point "<<i<<" initially was in "<<z_tree2[i]->name<<" and now is in "<<z_tree2_updated[i]->name<<endl;
}

cout<<"Get a vector that describes Tree 1"<<endl;
vector<int> desc1=G1->describeTree();
vec desc1_arma=stdToArmaVec(desc1); //if want arma
desc1_arma.print("Result:");

cout<<"Recreate the same tree using the created vector"<<endl;
Node* G1_copy=G1->translateTree(desc1); //can use any node to call translate tree

cout<<"Create a new tree of depth 2 and binary splits"<<endl;
Node* CreatedTree=G1->createTree(2,2);

cout<<"Create tree to practise swapping gates on"<<endl;
Gate* G_1=new Gate();
Gate* G_2=new Gate();
Gate* G_3=new Gate();
Expert* E_1=new Expert();
Expert* E_2=new Expert();
Expert* E_3=new Expert();
Expert* E_4=new Expert();
G_1->name="G_1";
G_2->name="G_2";
G_3->name="G_3";
E_1->name="E_1";
E_2->name="E_2";
E_3->name="E_3";
E_4->name="E_4";

G_1->addChild(E_1);
G_1->addChild(G_2);
G_2->addChild(E_2);
G_2->addChild(G_3);
G_3->addChild(E_3);
G_3->addChild(E_4);

G_1->Parent=NULL;
G_1->issueID();
G_1->issueIDLR();

E_1->expertmodel=NF;
E_2->expertmodel=NF;
E_3->expertmodel=NF;
E_4->expertmodel=NF;

vector<Node*> z_assign=assignPoints(G_1,10);

E_1->beta=E_1->expertmodel->findBeta(E_1->subsetY(y,E_1->getPointIndices(z_assign)),E_1->subsetX(X,E_1->getPointIndices(z_assign)),logsigma_sq,mu_beta,Sigma_beta);
E_2->beta=E_2->expertmodel->findBeta(E_2->subsetY(y,E_2->getPointIndices(z_assign)),E_2->subsetX(X,E_2->getPointIndices(z_assign)),logsigma_sq,mu_beta,Sigma_beta);
E_3->beta=E_3->expertmodel->findBeta(E_3->subsetY(y,E_3->getPointIndices(z_assign)),E_3->subsetX(X,E_3->getPointIndices(z_assign)),logsigma_sq,mu_beta,Sigma_beta);
E_4->beta=E_4->expertmodel->findBeta(E_4->subsetY(y,E_4->getPointIndices(z_assign)),E_4->subsetX(X,E_4->getPointIndices(z_assign)),logsigma_sq,mu_beta,Sigma_beta);
E_1->logsigma_sq=1;
E_2->logsigma_sq=1;
E_3->logsigma_sq=1;
E_4->logsigma_sq=1;

mat z_G_1=G_1->getZ(z_assign);
mat z_G_2=G_2->getZ(z_assign);
mat z_G_3=G_3->getZ(z_assign);

G_1->gamma=G_1->findGammaMLE(G_1->subsetX(X,G_1->getPointIndices(z_assign)),z_G_1,Omega);
G_2->gamma=G_2->findGammaMLE(G_2->subsetX(X,G_2->getPointIndices(z_assign)),z_G_2,Omega);
G_3->gamma=G_3->findGammaMLE(G_3->subsetX(X,G_3->getPointIndices(z_assign)),z_G_3,Omega);

vector<Node*> z_assign_new=G_1->MCMC(100,y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_assign);
for(int i=0;i<z_assign.size();i++){
    cout<<"Point "<<i<<" initially was in "<<z_tree1[i]->name<<" and now is in "<<z_assign_new[i]->name<<endl;
}

//G_1->updateSwap(G_2,G_3,1,y,X,z_assign);
//G_1->printDescendants();
//cout<<"The parent of "<<G_1->Children[1]->Children[0]->name<<endl;
//cout<<"is "<<G_1->Children[1]->Children[0]->Parent->name<<endl;


G_1->updateSwap(G_1,G_2,1,y,X,z_assign);
cout<<"descendants:"<<endl;
G_1->printDescendants();
//vec aa=stdToArmaVec(Trial->describeTree());
//aa.print("dscription:");
cout<<"descendants:"<<endl;
G_2->printDescendants();
delete G_2;
cout<<"descendants after deleting:"<<endl;
G_2->printDescendants();
cout<<"name after deleting:"<<endl;
cout<<G_2->name<<endl;



return 0;
}