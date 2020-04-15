//
// Created by 2022175V on 14/08/2019.
//
// Ludger is editing this file. Can you see this?

#include "Expert.h"
#include <iostream>
#include <vector>

using namespace std;
//using namespace arma;


/**
 * @brief Construct a new Expert:: Expert object.
 * 
 * @param aName a string name of the expert.
 */
Expert::Expert(string aName)
        :Node()
{
    this->type="E";
    this->name=aName;
    cout<<"Expert "<<name<<" has been created."<<endl;
}

/**
 * @brief Returns zero for the number of children of the expert.
 * 
 * @return int returns zero.
 */
int Expert::countChildren() {
    int n=0;
    return n;
}

/**
 * @brief Prints a message warning that experts can't have children.
 * 
 */
void Expert::printChildren() {
    cout<<"I am an expert "<<name<<", I can not have any children."<<endl;
}

/**
 * @brief Returns itself as a descendant.
 * 
 * @param desc vector of descendants of the expert.
 * @return vector<Node*> itself.
 */
vector<Node*> Expert::getDescendantsInternal(vector<Node*>* desc) {
    desc->push_back(this);
    return *desc;
}

/**
 * @brief Returns itself as a terminal node.
 * 
 * @param terminal vector of terminal nodes of the expert.
 * @return vector<Node*> itself.
 */
vector<Node*> Expert::getTerminalNodesInternal(vector<Node*>* terminal){
    terminal->push_back(this);
    return *terminal;
}

/**
 * @brief An internal function which helps to construct a numeric vector describing a tree. Returns zero for 
 * an expert.
 * This function is called at the node level and returns a vector containing the number of children 
 * of the gate. Then, the same is performed for each of the children of the gate. If the child is a gate,
 * the describeTreeInternal function will be called from the gate level. 
 * @param description vector to be filled with integers describing a tree.
 * @return vector<int> pointer to a vector with an entry of zero.
 */
vector<int> Expert::describeTreeInternal(vector<int>* description){
    description->push_back(0);
    return *description;
}
