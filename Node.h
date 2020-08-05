#ifndef MOE_NODE_H
#define MOE_NODE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

class Gate;

class Node{
    public:
    string name;
    int id;
    int idLR;
    vector<Node*>  Children;
    Gate* Parent;
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

};

#endif //MOE_NODE_H