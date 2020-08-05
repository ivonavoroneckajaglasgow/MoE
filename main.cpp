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
G3->addChild(E5);
G3->addChild(E6);

G1->issueIDLR();

int n=10;
int p=2;

mat X(n,p,fill::randn);
X.print("X:");

vector<Node*> all_experts=G1->getTerminalNodes();
vector<Node*> z_final(n);

for(int i=0;i<n;i++){
  int sub=rand() % all_experts.size();
  z_final[i]=all_experts[sub];
}

 cout<<"Allocations:"<<endl;

for(int i=0;i<z_final.size();i++){
   cout<<z_final[i]->name<<endl;
}

mat z_G1=G1->getZ(z_final);
z_G1.print("z for G1:");

mat z_G2=G2->getZ(z_final);
z_G2.print("z for G2:");

mat z_G3=G3->getZ(z_final);
z_G3.print("z for G3:");

vec G3_points=G3->getPointIndices(z_final);
G3_points.print("Points that went through G3:");

vec E4_points=E4->getPointIndices(z_final);
E4_points.print("Points that have been assigned to E4:");

vec E3_points=E3->getPointIndices(z_final);
E3_points.print("Points that have been assigned to E3:");


return 0; 
}