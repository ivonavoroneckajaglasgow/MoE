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
G2->gamma=G1->findGammaMLE(G2->subsetX(X,G2->getPointIndices(z_tree1)),z_G2,Omega);
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

G_1->issueID();
G_1->issueIDLR();

//Test copy constructor and delete children/parent:
Gate* G_1copy=G_1->copyThis();
Expert* E_1copy=E_1->copyThis();


vector<int> description=G_1->describeTree();
vec description_arma=stdToArmaVec(description);
description_arma.print("Describing vector for original G_1:"); //this will stay the same, just with a swap

vector<Node*> all_nodes=G_1->getDescendants();
for(int i=0;i<all_nodes.size();i++)
cout<<all_nodes[i]->name<<" "<<all_nodes[i]->idLR<<endl;
vec neworder("0 3 4 1 2 5");
vector<Node*> adj_nodes(neworder.size());
for(int i=0; i<neworder.size();i++){
adj_nodes[i]=all_nodes[neworder[i]];
cout<<adj_nodes[i]->name<<endl;
}

double loglik=G1->TotalLogLikelihood(y,X);
cout<<"Loglik of the whole tree "<<loglik<<endl;

vec somebeta2("1 1");
E_1->beta=somebeta2;


Gate* Copy=G_1->copyStructure();

vec somebeta("5 5");
cout<<"Beta of the expert child of copy"<<dynamic_cast<Expert*>(Copy->Children[0])->beta<<endl;
cout<<"Set it to be something else."<<endl;
dynamic_cast<Expert*>(Copy->Children[0])->beta=somebeta;
cout<<"Current beta of the expert child of copy"<<dynamic_cast<Expert*>(Copy->Children[0])->beta<<endl;
cout<<"Check if that affects E1"<<endl;
cout<<E_1->beta<<endl;

Gate* Swap=G_2->swap(G_3,1);


//We want a function that takes in the root node and returns:
//G_1->addChild(E_1);
//G_1->addChild(G_3);
//G_3->addChild(E_3);
//G_3->addChild(G_2);
//G_2->addChild(E_2);
//G_2->addChild(E_4);






return 0;
}