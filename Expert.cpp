#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "Expert.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Expert:: Expert object
 * 
 */
Expert::Expert(){
    cout<<"Expert has been created."<<endl;
}

/**
 * @brief Calculates Xbeta
 * 
 * @param X design matrix
 * @param beta parameter beta
 * @return vec result
 */
vec Expert::etafun(mat X, vec beta){
return X*beta;
}


/**
 * @brief Returns zero for the number of children of the expert
 * 
 * @return int zero
 */
int Expert::countChildren() {
    int n=0;
    return n;
}

/**
 * @brief Returns itself as a descendant
 * 
 * @param desc vector of descendants of the expert
 * @return vector<Node*> itself
 */
vector<Node*> Expert::getDescendantsInternal(vector<Node*>* desc) {
    desc->push_back(this);
    return *desc;
}

/**
 * @brief Returns itself as a terminal node
 * 
 * @param terminal vector of terminal nodes of the expert
 * @return vector<Node*> itself
 */
vector<Node*> Expert::getTerminalNodesInternal(vector<Node*>* terminal){
    terminal->push_back(this);
    return *terminal;
}

/**
 * @brief Sets the left to right ID 
 * 
 * @param start starting ID number
 * @return int incremented ID number
 */
int Expert::issueIDLR(int start){
    this->idLR=start;
    return start+1;
}

/**
 * @brief Returns its left to right ID 
 * Gate calls this function to get the right most ID
 * @return int left to right ID
 */
 int Expert::rightMostNodeID(){
     return this->idLR;
 }

/**
 * @brief Checks if node has the same left to right ID as this expert
 * 
 * @param node node to be checked
 * @return int yes/no 
 */
int Expert::isInRange(Node* node){
    if(this->idLR==node->idLR) return 1;
    else return 0;
}

/**
 * @brief Returns a matrix of ones
 * Used in path probability calculation, where multiplication by one does not change the result
 * @param X design matrix
 * @param gamma vector of gating parameters
 * @return mat matrix of ones
 */
mat Expert::pi_calculator(mat X, vec gamma){
    mat result(X.n_rows,gamma.size()/X.n_cols);
    result.fill(1);
    return result;
}

/**
 * @brief Updates parameters of the expert
 * 
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance value on log scale
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 */
void Expert::MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final){
    cout<<"Updating beta and sigma for Expert " <<this->name<<endl;
    //cout<<"Before: "<<this->beta<<endl;
    vec points=this->getPointIndices(z_final);
    //points.print("points:");
    vec myY=this->subsetY(y,points);
    mat myX=this->subsetX(X,points);
    this->beta=this->expertmodel->updateBeta(this->beta,myY,myX,logsigma_sq,mu_beta,Sigma_beta);
    this->logsigma_sq=this->expertmodel->updateSigma(this->logsigma_sq,myY,myX,this->beta,a,b,static_cast<int>(points.size()));
    //cout<<"After: "<<this->beta<<endl;
}

/**
 * @brief Function that tanslates current state of the tree to a string in json format
 * 
 * @param indent spacing variable 
 * @return string describing current state of the expert
 */
 string Expert::jsonify(int indent) {
  map<string, string> m;
  map<string, string> p;
  string s;
  m["__name__"] = "\""+this->name+"\"";
  m["__type__"] = "\"expert\"";
  m["__family__"] = "\"normal\"";
  m["_beta_"] = "\n" + vec2arraystring(this->beta, indent+6);
  ostringstream os;
  os << scientific << this->logsigma_sq;
  p["_logsigma_sq_"] = os.str();
  m["_par_"] = "\n" + jsondict(p, indent+6);
  return jsondict(m, indent);
} 

