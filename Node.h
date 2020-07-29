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
    vector<Node*>  Children;
    Gate* Parent;
    void printParent();
    virtual void printDescendants();
    virtual void printTerminalNodes();
};

#endif //MOE_NODE_H