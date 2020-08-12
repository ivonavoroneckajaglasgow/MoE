#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Expert.h"

using namespace std;
using namespace arma;

Expert::Expert(){
    cout<<"Expert has been created."<<endl;
}

vec Expert::etafun(mat X, vec beta){
return X*beta;
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

int Expert::issueIDLR(int start){
    this->idLR=start;
    return start+1;
   // return -99;
}

 int Expert::rightMostNodeID(){
     return this->idLR;
 }

int Expert::isInRange(Node* node){
    if(this->idLR==node->idLR) return 1;
    else return 0;
}

mat Expert::pi_calculator(mat X, vec gamma){
    mat result(X.n_rows,gamma.size()/X.n_cols);
    result.fill(1);
    return result;
}