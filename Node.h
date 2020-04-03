//
// Created by 2022175V on 14/08/2019.
//

#ifndef MOE_3_NODE_H
#define MOE_3_NODE_H

#include <iostream>
#include <vector>

using namespace std;

class Gate;
class Expert;

class Node {
public:
    int  id;
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
};


#endif //MOE_3_NODE_H