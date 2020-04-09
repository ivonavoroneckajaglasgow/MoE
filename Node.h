//
// Created by 2022175V on 14/08/2019.
//

#ifndef MOE_3_NODE_H
#define MOE_3_NODE_H

#include <iostream>
#include <vector>

#include "GateParameters.h"
#include "NormalParameters.h"

using namespace std;
using namespace arma;

class Gate;
class Expert;

class Node {
public:
    int  id;
    Node* helper;
    Gate* Parent;
    string name;
    string type;
    Node();
    void printParent();
    virtual void printChildren();
    void printAncestors();

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
    vector<int> describeTree();
    virtual vector<int> describeTreeInternal(vector<int>* description);
    virtual void issueID();
    virtual void issueID_helper1(int* gate_id, int* expert_id);
    virtual void issueID_helper2(int* gate_id, int* expert_id);
    static Node* createTreeInternal(GateParameters aParameters,NormalParameters aParameters2, int depth, int nchildren, int* gcount, int* ecount);
    Node* createTree(GateParameters aParameters,NormalParameters aParameters2, int depth, int nchildren);
    int populateGate(Gate* parent, vector<int> description, int start, int* gcount, int* ecount, GateParameters aParameters, NormalParameters aParameters2);
    Node* translateTree(vector<int> description, GateParameters aParameters, NormalParameters aParameters2);
    };


#endif //MOE_3_NODE_H