#ifndef MOE_NODE_H
#define MOE_NODE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

class Gate;
class ExpertModel;

class Node{
    public:
    string name; //name 
    int id;      //top to bottom ID
    int idLR;    //left to right ID 
    vector<Node*>  Children;//vector of pointers to children
    Gate* Parent;           //pointer to parent
    virtual Node* copyThis();         //creates a deep copy of the object
    void printParent();               //prints the name of the parent
    virtual void printChildren();
    virtual void printDescendants();  //prints names of descendants 
    virtual void printTerminalNodes();//prints names of terminal nodes
    virtual Gate* getParent();//returns pointer to the parent node
    virtual vector<Node*> getChildren();//returns vector of pointers to children
    vector<Node*> getAncestors();//returns vector of pointers to ancestors
    vector<Node*> getAncestorsInternal(vector<Node*>* ancest);//helper for the above function
    vector<Node*> getDescendants();//returns vector of pointers to descendants 
    virtual vector<Node*> getDescendantsInternal(vector<Node*>* desc);//helper for the above function
    vector<Node*> getTerminalNodes();//returns a vector of pointers to the terminal nodes 
    virtual vector<Node*> getTerminalNodesInternal(vector<Node*>* terminal);//helper for the above function
    virtual int countChildren();//returns the total number of children
    virtual void issueID();//issues top to bottom ID
    virtual void issueID_helper1(int* gate_id, int* expert_id);//helper for the above
    virtual void issueID_helper2(int* gate_id, int* expert_id);//helper for the above
    virtual int issueIDLR(int start);//issues left to right ID
    int leftMostNodeID();//returns left most node ID
    virtual int rightMostNodeID();//returns right most node ID
    virtual int isInRange(Node* node);//checks if node is in the range of descendants
    vec getDescendantRange(); //returns a vector of descendant range left to right IDs
    int countPoints(vector<Node*> z_final);//counts how many points have travelled through node
    vec getPointIndices(vector<Node*> z_final);//return a vector of indeces for points that have travelled through node
    mat subsetX(mat X, vec index);//subsets a matrix X
    vec subsetY(vec Y, vec index);//subsets a vector y
    virtual mat pi_calculator(mat X, vec gamma);//calculates mixing proportions
    virtual void MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final);//updates parameters of a node 
    virtual mat getZ(vector<Node*> z_final);//returns a matrix of allocations z
    virtual string jsonify(int indent);//writes the current state of the tree as a json format string
    vector<int> describeTree(); //describes the tree as a vector of integers
    virtual vector<int> describeTreeInternal(vector<int>* description); //helper for the above
    int populateGate(Gate* parent, vector<int> description, int start, int* gcount, int* ecount); //populates the supplied gate with children
    Node* translateTree(vector<int> description); //translates a vector of integers into a tree
    static Node* createTreeInternal(int depth, int nchildren, int* gcount, int* ecount); //wrapper for the function below
    Node* createTree(int depth, int nchildren); //creates a tree object given a set of instructions
    Gate* mostSeniorGate();
    Node* findNode(string Name);
};

string jsondict(map<string, string> m, int indent);
string vec2arraystring(vec b, int indent);
string mat2arraystring(Mat<double> A, int indent);

#endif //MOE_NODE_H