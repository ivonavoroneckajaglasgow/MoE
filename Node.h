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
    string name;
    int id;
    int idLR;
    vector<Node*>  Children;
    Gate* Parent;
    //vec beta; //does nothing only here because terminal nodes are Node*
    //double logsigma_sq;   //does nothing only here because terminal nodes are Node* 
    //ExpertModel* expertmodel; //does nothing only here because terminal nodes are Node*
    void printParent();
    virtual void printDescendants();
    virtual void printTerminalNodes();
    virtual Gate* getParent();
    virtual vector<Node*> getChildren();
    vector<Node*> getAncestors();
    vector<Node*> getAncestorsInternal(vector<Node*>* ancest);
    vector<Node*> getDescendants();
    virtual vector<Node*> getDescendantsInternal(vector<Node*>* desc);
    vector<Node*> getTerminalNodes();
    virtual vector<Node*> getTerminalNodesInternal(vector<Node*>* terminal);
    virtual int countChildren();
    virtual void issueID();
    virtual void issueID_helper1(int* gate_id, int* expert_id);
    virtual void issueID_helper2(int* gate_id, int* expert_id);
    virtual int issueIDLR(int start);
    int leftMostNodeID();
    virtual int rightMostNodeID();
    virtual int isInRange(Node* node);
    vec getDescendantRange(); 
    int countPoints(vector<Node*> node);
    vec getPointIndices(vector<Node*> node);
    mat subsetX(mat X, vec index);
    vec subsetY(vec Y, vec index);
    //virtual mat getPathProbs(mat X, vector<Node*> z_final);
    virtual mat pi_calculator(mat X, vec gamma);
    virtual void MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final);
    virtual mat getZ(vector<Node*> z_final);

};

#endif //MOE_NODE_H