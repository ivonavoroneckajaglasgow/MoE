#ifndef MOE_GATE_H
#define MOE_GATE_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"
#include "Expert.h"
#include "Node.h"

using namespace std;
using namespace arma;


class Gate: public Node{
    public:
    Gate();    //constructor
    Gate(const Gate &gate); //copy constructor
    vec gamma; //parameter associated with this gate length rp
    mat Omega; //prior variance for gamma
    Gate* copyThis();            //creates a new gate that is a copy of this one
    void copyStructure(Gate* Copy);//creates a copy of the gate and all its descedants and stores that in Copy
    Gate* copyStructure();       //wrapper for the above
    void addChild(Node* aChild); //function that adds a child to the gate and gate as a parent to the child
    void deleteChildren();       //deletes all children and itself as a parent
    void printChildren();        //prints children names 
    void printDescendants();     //prints descendant names
    void printTerminalNodes();   //prints termial nodes (experts) names
    vector<Node*> getChildren(); //returns a vector with pointers to children
    vector<Node*> getDescendantsInternal(vector<Node*>* desc);       //returns a vector with pointers to descendants
    vector<Node*> getTerminalNodesInternal(vector<Node*>* terminal); //returns a vector with pointers to termial nodes (experts)
    int countChildren();   //returns the total number of children
    int countDescendants();//returns the total number of descendants
    void issueID();//issues IDs top to bottom with the help of two functions below
    void issueID_helper1(int* gate_id, int* expert_id);//helper for issueID
    void issueID_helper2(int* gate_id, int* expert_id);//helper for issueID
    double loglik(mat z, mat pi);//the log-likelihood for the model given allocations z and mixing proportions pi
    double loglik(mat X, vec gamma, mat z);//the log-likelihood for the model given the design matrix X, gating parameters gamma and allocations z
    mat pi_calculator(mat X, vec gamma); //pi calculator that saves intermediate results
    mat pi_calculator2(mat X, vec gamma);//pi calculator that loops over each observation
    mat pi_calculator3(mat X, vec gamma);//pi calculator that subtracts the maximum value
    vec score(mat X, mat z, mat pi); //score for IWLS
    mat hessian(mat X, mat pi);      //hessian for IWLS
    vec findGammaMLEChol(mat X, mat z, mat Omega, mat* L); //estimate gamma using: Cholesky decomposition + get poiter to L matrix
    vec findGammaMLEChol(mat X, mat z, mat Omega);         //estimate gamma using: Cholesky decomposition
    vec findGammaMLEQR(mat X, mat z, mat Omega, mat* R);   //estimate gamma using: QR decomposition + get pointer to R matrix
    vec findGammaMLEQR(mat X, mat z, mat Omega);           //estimate gamma using: QR decomposition
    vec findGammaMLE(mat X, mat z, mat Omega);             //estimate gamma using: no decomposition
    mat makeAchol(vec pi);     //performs Cholesky decomposition on a matrix required for QR decomposition in findGammaQR() below
    mat getXout(mat X, mat pi);//outputs big matrix X required in findGammaQR() below
    vec getZeta(mat z, mat pi);//outputs the new response required in findGammaQR() below
    vec findGammaQR(mat X, mat z, mat Omega, mat* R); //estimates gamma using adapted IWLS and QR decomposition + get pointer to R matrix in the QR decomposition
    vec findGammaQR(mat X, mat z, mat Omega);         //estimates gamma using adapted IWLS and QR decomposition
    vec updateGamma(vec gammaold, mat X, mat z, mat Omega); //proposes a new value for gamma and checks if it is accepted
    vec logmvndensity(vec response, vec mean, mat Sigma);   //multivariate normal density function on a log-scale
    vec logmvndensity(vec response, vec mean, mat* R);      //multivariate normal density function on a log-scale + takes a pointer to R from QR decomp as an argument instead of Sigma
    mat getRowSumsMat(mat A);  //calculates rowSums
    int issueIDLR(int start);  //issues ID left to right
    int issueIDLR();           //wrapper for the above function
    int rightMostNodeID();     //returns the right most node ID
    int isInRange(Node* node); //checks if node is in the range of descendant IDs
    rowvec getZ_perpoint(Node* node); //returns a vector of length of number of splits, entry of 1 indicates which split node is in
    mat getZ(vector<Node*> z_final); //returns allocations matrix z
    vector<Node*> updateZ(vec y, mat X, vector<Node*> z_final); //returns a vector of updated allocations as pointers to experts 
    vec getSampleProbsForZ(vec y, mat X); //X and y are for one i only, returns alpha vector containing of log(pi)+logdensity
    Node* updateZ_onepoint_sample(vec alpha);//allocates a point to an expert with supplied probabilities
    double getPathProb(Node* node, mat X); //X is for one i, returns path probability down from this to node
    double getPathProb_internal(Node* node, mat X, double result); //helper for the above function
    vec getPathProb_mat(Node* node, mat X); //returns path probabilities (rows are observations and columns are experts)
    vec getChildrenIndicesLR(); //returns LR indices of children
    int whichChild(Node* node); //returns which child node is (to help identify the reference split)
    vec predict(mat X);//returns y hat
    void MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final); //updates gamma for a gate 
    vector<Node*> MCMC_OneRun(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final);//updates parameters followed by allocations once 
    vector<Node*> MCMC(int N, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final);//updates parameters followed by allocations N times 
    string jsonify(int indent); //produces a JSON string describing the current state
    string jsonify();           //wrapper for the above function
    vector<Node*> updateZparamsIntegratedOut(vec y, mat X, vector<Node*> z_final,vec mu_beta, mat Sigma_beta, double a, double b);//updates allocations z with parameters beta and sigma integrated out (Normal Model only)
    vector<int> describeTreeInternal(vector<int>* description);//helps describe the tree as a vector of integers
    double TotalLogLikelihood(vec y, mat X); //calculates the brute force log likelihood for the whole tree
    Gate* swap(Gate* gate, int replace); //swaps the replace'th child of gate by this
    };
    

#endif //MOE_GATE_H