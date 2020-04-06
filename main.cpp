#include <iostream>
#include <vector>
#include "Gate.h"
#include "Expert.h"
#include "NormalExpert.h"
#include "NormalParameters.h"
#include "GateParameters.h"

using namespace std;

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

/**
 * @brief Creates a tree object given a set of instructions
 * Creates a node, which can be either an expert (if depth is zero) or a gate
 * @param depth sets the depth of the tree to be created
 * @param nchildren sets the number of children to add at each split of the tree
 * @param gcount pointer to a variable, which counts/tracks the number of gates in the tree
 * @param ecount pointer to a variable, which counts/tracks the number of experts in the tree
 * @return Node* pointer to the root node of the newly created tree
 */
Node* createTree(int depth, int nchildren, int* gcount, int* ecount) {    
    if (depth==0)
        return new NormalExpert("E" + std::to_string((*ecount)++), generalNormalParams);
    Gate* root = new Gate("G" + std::to_string((*gcount)++), generalGateParams);
    for (int i=0; i<nchildren; i++)
        root->addChild(createTree(depth-1, nchildren, gcount, ecount));
    return root;
}

/**
 * @brief Internal function, used in the translateTree function
 * Populates the supplied gate with children. Decides on whether to add a Gate or an Expert as a child based
 * on the supplied description vector
 * @param parent the gate to be populated
 * @param description vector of integers containing the number of children of each node  
 * @param start integer that denotes the position of the description vector
 * @param gcount pointer to a variable, which counts/tracks the number of gates in the tree 
 * @param ecount pointer to a variable, which counts/tracks the number of experts in the tree 
 * @return int returns the integer denoting last position in the description vector
 */
int populateGate(Gate* parent, vector<int> description, int start, int* gcount, int* ecount) {
    int pos =  start;
    for (int i=0; i<description[start]; i++) {
        pos++;
        if (pos >= description.size() || description[pos] == 0) {
            parent->addChild(new NormalExpert("E" + std::to_string((*ecount)++), generalNormalParams));
        } else {
            Gate *g = new Gate("G" + std::to_string((*gcount)++), generalGateParams);
            parent->addChild(g);
            pos = populateGate(g, description, pos, gcount, ecount);
        }
    }
    return pos;
}

/**
 * @brief Turns a numerical vector describing a tree into a tree
 * Turns a numerical vector of integers (as produced by Gate->describeTree()) into a tree object
 * @param description vector of integers containing the number of children of each node
 * @return Node* pointer to the root node of the newly created tree
 */
Node* translateTree(vector<int> description) {
    if(accumulate(description.begin(), description.end(), 0)!=description.size()-1)
        cout<<"Warning: the description vector is not complete. All missing entries will be replaced by 1, i.e. an Expert will be added."<<endl;
      // throw "Warning: the description vector is not complete. All missing entries will be replaced by 1, i.e. an Expert will be added.";
    int ecount=0;
    int gcount=0;
    if (description.size()==0)
        return 0;
    if (description[0]==0)
        return new NormalExpert("E" + std::to_string(ecount++), generalNormalParams);
    Gate* root = new Gate("G" + std::to_string(gcount++), generalGateParams);
    populateGate(root, description, 0, &gcount, &ecount);
    return root;
}

int main(){

/* Gate* G = new Gate("G",generalGateParams);
NormalExpert* E1 = new NormalExpert("E1",generalNormalParams);
NormalExpert* E2 = new NormalExpert("E2",generalNormalParams);

G->addChild(E1);
G->addChild(E2);

vector<Node*> test;
test=G->getDescendants();
for(int i=0;i<test.size();i++)
cout<<test[i]->name<<endl; */

     int  n_g=0;
     int  n_e=0;
     int  d=2;
     int  chl=2;

/* Gate* Gtest= new Gate("Gtest",generalGateParams, generalNormalParams,d,chl,&n_g,&n_e);

cout<<Gtest->name<<endl;
Gtest->printChildren(); */

//      int n_g2=0;
//      int n_e2=0;
//      int  d2=2;
//      int  chl2=2;

// Gtest->createTree(generalGateParams, generalNormalParams,d2,chl2,&n_g2,&n_e2);

return 0;
}

