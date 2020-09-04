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
    vec* y;       //data observed
    vec  eta;     // XB
    mat* X;
    vec beta;
    double logsigma_sq;  // some variance parameter  
    ExpertModel* expertmodel; // the model for this expert
    Expert();//constructor
    vec etafun(mat X, vec beta);
    int countChildren();
    vector<Node*> getDescendantsInternal(vector<Node*>* desc);
    vector<Node*> getTerminalNodesInternal(vector<Node*>* terminal);
    int issueIDLR(int start);
    int rightMostNodeID();
    int isInRange(Node* node);
    mat pi_calculator(mat X, vec gamma);
    void MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final);
    string createJSON();
    string createJSON2 (string s);
    string jsonify(int indent);
};

#endif //MOE_EXPERT_H