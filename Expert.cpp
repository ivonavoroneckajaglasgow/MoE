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

void Expert::MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final){
    cout<<"Updating beta and sigma for Expert " <<this->name<<endl;
    //cout<<"Before: "<<this->beta<<endl;
    vec points=this->getPointIndices(z_final);
    vec myY=this->subsetY(y,points);
    mat myX=this->subsetX(X,points);
    this->beta=this->expertmodel->updateBeta(this->beta,myY,myX,logsigma_sq,mu_beta,Sigma_beta);
    this->logsigma_sq=this->expertmodel->updateSigma(this->logsigma_sq,myY,myX,this->beta,a,b,points.size());
    //cout<<"After: "<<this->beta<<endl;
}

string Expert::createJSON(){
    
//Creating JSON output for expert//
//create a template string
    string ExpertJSON=("{type: Expert beta: [BTA] family: FML additional_parameters: ADPAR}");
//change the FML part of the string to the family of the expert in question. Will need to create 
// a string method that outputs "Gaussian","Binomial" etc based on the expertmodel.
    ExpertJSON.replace(ExpertJSON.find("FML"),3,"Gaussian");

//Turn beta into a string
    vec myvec=this->beta;  
    ostringstream ss;

// .st() to transpose column vector into row vector
    myvec.st().raw_print(ss);  

// get string version of vector with an end-of-line character at end
    string s1 = ss.str();

// remove the end-of-line character at end
    string mystring = s1.substr(0, (s1.size() > 0) ? (s1.size()-1) : 0);

//replace BTA by the string containing beta */
    ExpertJSON.replace(ExpertJSON.find("BTA"),3,mystring);

//finally, if the family is Gaussian replace ADPAR by additional parameters
    ostringstream s_2;
    s_2 <<  this->logsigma_sq;
    string mystring_2 = s_2.str();

    ExpertJSON.replace(ExpertJSON.find("ADPAR"),5,mystring_2);

    return ExpertJSON;

}