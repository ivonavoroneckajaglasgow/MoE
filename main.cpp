#include <iostream>
#include <vector>
#include "Gate.h"
#include "Expert.h"
#include "NormalExpert.h"
#include "NormalParameters.h"
#include "GateParameters.h"

using namespace std;
using namespace arma;

NormalParameters generalNormalParams{
        {0.0,0.0},
        1.0,
        {0.0,0.0},
        {{1.0,0.0},{0.0,2.0}},
        0.001,
        0.001
};

GateParameters generalGateParams{
        {0.0,0.0},
        {3.0,10.0},
        {{10.0,0.0},{0.0,10.0}}
};

int main(){

Node* my_node;
int d=2;
int chl=2;
Node* new_node;
new_node=my_node->createTree(generalGateParams,generalNormalParams, d,chl);
vector<int> desc;
desc=new_node->describeTree();
for(int i=0;i<desc.size();i++)
cout<<desc[i]<<endl;

Node* new_node_recreated;
new_node_recreated=my_node->translateTree(desc,generalGateParams,generalNormalParams);

return 0;
}

