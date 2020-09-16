#ifndef MOE_EXPERT_H
#define MOE_EXPERT_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "ExpertModel.h"
#include "Node.h"

using namespace std;
using namespace arma;

//This is a superclass of NormalExert and ExpertModel objects 
//Most functions are virtual and overwritten at subclass levels

class Expert: public Node{
    public:
    vec* y;            //pointer to data in the Data object
    mat* X;            //pointer to data in the Data object
    vec  eta;          //XB
    vec beta;          //beta parameter
    double logsigma_sq;//variance parameter  
    ExpertModel* expertmodel; //the model for this expert
    Expert();                   //constructor
    vec etafun(mat X, vec beta);//eta calculator
    int countChildren(); //counts the number of children, which is zero for experts
    vector<Node*> getDescendantsInternal(vector<Node*>* desc); //returns itself as a descendant 
    vector<Node*> getTerminalNodesInternal(vector<Node*>* terminal);//returns itself as a terminal node
    int issueIDLR(int start); //issues ID left to right
    int rightMostNodeID();  //returns its own ID and Gate decides which one is right most
    int isInRange(Node* node);//returns yes or no to whether its ID is te same as node
    mat pi_calculator(mat X, vec gamma);//returns a matrix of ones (used in multiplication when calculating path probabilities)
    void MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final);//updates beta and sigma
    string jsonify(int indent); //produces a JSON string describing the current state
    vector<int> describeTreeInternal(vector<int>* description);//helps describe the tree as a vector of integers
};

#endif //MOE_EXPERT_H