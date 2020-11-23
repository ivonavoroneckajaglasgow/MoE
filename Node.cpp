#define _USE_MATH_DEFINES
#define EPS 2.22e-16

#include <vector>
#include <cmath>

#include "armadillo"
#include <boost/algorithm/string/join.hpp>

using namespace std;
using namespace arma;

#include "Node.h"
#include "Gate.h"

/**
 * @brief Prints the name of the parent
 * 
 */
void Node::printParent() {
        cout << name << " parent is " << Parent->name << "." << endl;
}

/**
 * @brief Prints names of descendants
 * Handled at Gate and Expert levels
 */
void Node::printDescendants(){

}

/**
 * @brief Prints names of terminal nodes
 * Handled at Gate and Expert levels
 */
void Node::printTerminalNodes(){

}

/**
 * @brief Creates a deep copy of the object
 * Handled at Gate and Expert levels
 * @return Node* pointer to the object copy
 */
Node* Node::copyThis(){
    return 0;
}

/**
 * @brief Returns a pointer to parent node
 * 
 * @return Gate* pointer to parent node
 */
Gate* Node::getParent(){
    return Parent;
}

/**
 * @brief Returns a vector of pointers to children
 * Handled at Gate and Expert levels
 * @return vector<Node*> empty vector
 */
vector<Node*> Node::getChildren() {
    vector<Node*> empty;
    return empty;
}

/**
 * @brief Returns ancestors of the node
 * Wrapper for the internal function getAncestorsInternal()
 * @return vector<Node*> vector of pointers to the ancestors of the node
 */
vector<Node*> Node::getAncestors() {
    vector <Node*> ancest_test;
    ancest_test=this->getAncestorsInternal(&ancest_test);
    return ancest_test;
}

/**
 * @brief An internal helper function that retrieves ancestors of the node
 * 
 * @param ancest vector to be filled in with ancestors of the node
 * @return vector<Node*> vector of pointers to the ancestors of the node
 */
vector<Node*> Node::getAncestorsInternal(vector<Node*>* ancest) {

    if(this->Parent!=NULL){
        ancest->push_back(this->Parent);
        this->Parent->getAncestorsInternal(ancest);
    }
    return *ancest;

}

/**
 * @brief Returns the descendents of the node
 * Calls the getDescendantsInternal() function from gate and expert levels
 * @return vector<Node*> vector of pointers to the descendants of the node
 */
vector<Node*> Node::getDescendants() {
    vector <Node*> desc_test;
    desc_test=this->getDescendantsInternal(&desc_test);
    return desc_test;
}

/**
 * @brief Returns a vector of pointers to the descendants, wrapper for getDescendants()
 * Handled at expert and node levels
 * @param desc vector to be filled in with pointers to the descendants
 * @return vector<Node*> vector of pointers to the descendants
 */
vector<Node*> Node::getDescendantsInternal(vector<Node*>* desc){
vector<Node*> a;
return a;
}

/**
 * @brief Returns the vector of pointers to the terminal nodes, wrapper for getTerminalNodesInternal()
 * Calls the getTerminalNodesInternal() function from gate and expert levels
 * @return vector<Node*> vector of pointers to the terminal nodes 
 */
vector<Node*> Node::getTerminalNodes() {
    vector <Node*> terminal_test;
    terminal_test=this->getTerminalNodesInternal(&terminal_test);
    return terminal_test;
}

/**
 * @brief Returns the vector of pointers to the terminal nodes
 * Handled at expert and node levels
 * @param terminal vector to be filled in by pointers to the terminal nodes
 * @return vector<Node*> vector of pointers to the terminal nodes
 */
vector<Node*> Node::getTerminalNodesInternal(vector<Node*>* terminal){
vector<Node*> a;
return a;
}

/**
 * @brief Returns the total number of children
 * Handled at Gate and Expert levels
 * @return int number of children
 */
int Node::countChildren() {
return 0;
}

/**
 * @brief Issues top to bottom IDs
 * Handled at Gate and Expert levels
 */
void Node::issueID() {}

/**
 * @brief Issues top to bottom IDs
 * Handled at Gate and Expert levels
 * @param gate_id tracked of gate IDs
 * @param expert_id trackr of expert IDs
 */
void Node::issueID_helper1(int* gate_id, int* expert_id){}

/**
* @brief Issues top to bottom IDs
 * Handled at Gate and Expert levels
 * @param gate_id tracked of gate IDs
 * @param expert_id trackr of expert IDs 
 */
void Node::issueID_helper2(int* gate_id, int* expert_id) {}

/**
 * @brief Issues left to right IDs
 * Handled at Gate and Expert levels
 * @param start starting ID value
 * @return int fiishing ID value
 */
int Node::issueIDLR(int start){
    return 0;
}

/**
 * @brief Returns left most descendant ID
 * 
 * @return int left to right left most ID
 */
int Node::leftMostNodeID(){
    return this->idLR;
}

/**
 * @brief Returns right most descendant ID
 * Handled at Gate and Expert levels
 * @return int left to right right most ID
 */
int Node::rightMostNodeID(){
    return 0;
}

/**
 * @brief Checks if node is in the range of descendants
 * 
 * @param node node to check
 * @return int yes/no
 */
int Node::isInRange(Node* node){
    return 0;
}

/**
 * @brief Returns a vector containing the range of left to right IDs of descendants
 * 
 * @return vec vector containing the range of left to right IDs of descendants
 */
vec Node::getDescendantRange(){
    vec result;
    result<<this->leftMostNodeID()<<this->rightMostNodeID();
    return result;
}


/**
 * @brief Counts the number of points that have been assigned to the node based on final allocations
 * 
 * @param node a vector of pointers to experts (the final allocations)
 * @return int the number of points that have been assigned to the node
 */
int Node::countPoints(vector<Node*> z_final){
   int count=0;
   for(int i=0;i<z_final.size();i++){
    if(this->isInRange(z_final[i])==1){
       count=count+1;
    }
   }
   return count;
}

/**
 * @brief Outputs a vector containing indeces of the points that have travelled through the node based on final allocations
 * 
 * @param z_final a vector of pointers to experts (the final allocations)
 * @return vec vector containing indeces of the points that have travelled through the node
 */
vec Node::getPointIndices(vector<Node*> z_final){
    vec range=this->getDescendantRange();    
    vec result(this->countPoints(z_final));
    int position=0;
    for(int i=0;i<z_final.size();i++){
        if(this->isInRange(z_final[i])==1){
            result[position]=i;
            position=position+1;
        }
    }
    return result;
}

/**
 * @brief Subsets a matrix X
 * 
 * @param X matrix
 * @param index row indeces to be subset
 * @return mat resulting matrix
 */
mat Node::subsetX(mat X, vec index){
    mat result(index.size(),X.n_cols);
    for(int i=0;i<index.size();i++){
        result.row(i)=X.row(static_cast<int>(index[i]));
    }
    return result;
}

/**
 * @brief Subsets a vector y
 * 
 * @param y vector 
 * @param index point indeces to be subset
 * @return vec resulting vector
 */
vec Node::subsetY(vec y, vec index){
    vec result(index.size());
    for(int i=0;i<index.size();i++){
        result[i]=y[static_cast<int>(index[i])];
    }
    return result;
}

/**
 * @brief Returns a matrix of mixing proportions
 * Handled at Gate and Expert levels
 * @param X design matrix
 * @param gamma vector of gating parameters
 * @return mat matrix of mixing proportions
 */
mat Node::pi_calculator(mat X, vec gamma){
    return 0;
}

/**
 * @brief Updates parameters of the node
 * Handled on Gate and Expert levels
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq variance value on log scale
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 */
void Node::MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final){
    
}

/**
 * @brief returns allocations matrix z
 * Handled on Gate and Expert levels 
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 * @return mat allocations matrix z
 */
mat Node::getZ(vector<Node*> z_final){
    return 0;
}

/**
 * @brief Function that tanslates current state of the tree to a string in json format
 * Handled on Gate and Expert levels 
 * @param indent spacing variable (always 0 see wrapper below) 
 * @return string describing current state of the tree
 */
 string Node::jsonify(int indent){
    return string("?");
}

/**
 * @brief 
 * 
 * @param m 
 * @param indent 
 * @return string 
 */
string jsondict(map<string, string> m, int indent) {
    ostringstream os; 
    map<string, string>::iterator it;
    os << string(indent, ' ') << "{";    
    string s = string("");
    for (it = m.begin(); it!=m.end(); it++ ) {        
        os  << s << endl << string(indent+2, ' ') << "\"" << it->first << "\": " << it->second;
        s = string(",");
    }
    os << endl << string(indent, ' ') << "}";
    return os.str();
}

/**
 * @brief 
 * 
 * @param dv 
 * @return string 
 */
string dv2csv(vector<double> dv) {
  vector<string> sv;
  transform(begin(dv),
                 end(dv),
                 back_inserter(sv),
                 [](double d) {
                   ostringstream dos;
                   dos << scientific << d;
                   return dos.str();
                 });
  return boost::algorithm::join(sv, ", ");
}

/**
 * @brief 
 * 
 * @param b 
 * @param indent 
 * @return string 
 */
string vec2arraystring(vec b, int indent)
{
  vector<double> dv;
  ostringstream os; 
  dv = conv_to<vector<double>>::from(b);
  os << string(indent, ' ') << "[ " << dv2csv(dv) << " ]";
  return os.str();
}

/**
 * @brief 
 * 
 * @param A 
 * @param indent 
 * @return string 
 */
string mat2arraystring(Mat<double> A, int indent)
{
  vector<double> dv;
  ostringstream os;
  for (size_t i = 0; i < A.n_rows; ++i)
  {
    dv = conv_to<vector<double>>::from(A.row(i));
    if (i == 0)
    {
      os << string(indent, ' ') << "[ ";
    }
    else
    {
      os << string(indent + 2, ' ');
    }
    os << "[ " << dv2csv(dv) << " ]";
    if (i == A.n_rows - 1)
    {
      os << " ]";
    }
    else
    {
      os << "," << endl;
    }
  }
  return os.str();
} 

/**
 * @brief Function which constructs a numeric vector describing a tree
 * This function returns a vector containing the number of children of each node in the tree
 * To do so, the describeTreeInternal() is called on expert and gate level
 * @param description_test vector to be filled with integers describing a tree
 * @return vector<int> pointer to a vector of integers with the number of children of the gate
 */
vector<int> Node::describeTree(){
    vector <int> description_test;
    description_test=this->describeTreeInternal(&description_test);
    return description_test;
}

/**
 * @brief Function which helps construct a numeric vector describing a tree 
 * Handled at Expert and Node levels
 * @param description vector to be filled with integers describing a tree
 * @return vector<int> pointer to a vector of integers with the number of children of the gate
 */
vector<int> Node::describeTreeInternal(vector<int>* description){
   vector<int> a;
   return a;
}

/**
 * @brief Turns a numerical vector describing a tree into a tree
 * Turns a numerical vector of integers (as produced by Gate->describeTree()) into a tree object
 * @param description vector of integers containing the number of children of each node
 * @return Node* pointer to the root node of the newly created tree
 */
Node* Node::translateTree(vector<int> description) {
    //if(accumulate(description.begin(), description.end(), 0)!=description.size()-1)
      //  cout<<"Warning: the description vector is not complete. All missing entries will be replaced by 1, i.e. an Expert will be added."<<endl;
    int ecount=1;
    int gcount=1;
    if (description.size()==0)
        return 0;
    if (description[0]==0){
        Expert* e= new Expert();
        e->name="E" + std::to_string(ecount++);
        return e;
    }
    Gate* root = new Gate();
    root->name="G" + std::to_string(gcount++);
    this->populateGate(root, description, 0, &gcount, &ecount);
    return root;
}

/**
 * @brief Internal function, used in the translateTree() function
 * Populates the supplied gate with children. Decides on whether to add a Gate or an Expert as a child based
 * on the supplied description vector
 * @param parent the gate to be populated
 * @param description vector of integers containing the number of children of each node
 * @param start integer that denotes the position of the description vector
 * @param gcount pointer to a variable, which counts/tracks the number of gates in the tree
 * @param ecount pointer to a variable, which counts/tracks the number of experts in the tree
 * @return int returns the integer denoting last position in the description vector
 */
int Node::populateGate(Gate* parent, vector<int> description, int start, int* gcount, int* ecount) {
     int pos =  start;
    for (int i=0; i<description[start]; i++) {
        pos++;
        if (pos >= description.size() || description[pos] == 0) {
            Expert* e= new Expert();
            e->name="E" + std::to_string((*ecount)++);
            parent->addChild(e);
        } else {
            Gate* g = new Gate();
            g->name="G" + std::to_string((*gcount)++);
            parent->addChild(g);
            pos = populateGate(g, description, pos, gcount, ecount);
        }
    }
     return pos;
}

/**
 * @brief A wrapper for createTreeInternal
 * Initialises gcount and ecount automatically
 * @param depth sets the depth of the tree to be created
 * @param nchildren sets the number of children to add at each split of the tree
 * @return Node* pointer to the root node of the newly created tree
 */

Node* Node::createTree(int depth, int nchildren){
    int gcount=0;
    int ecount=0;
    Node* root;
    root=this->createTreeInternal(depth,nchildren,&gcount,&ecount);
    return root;
}
/**
 * @brief Creates a tree object given a set of instructions
 * Creates a node, which can be either an expert (if depth is zero) or a gate
 * @param depth sets the depth of the tree to be created
 * @param nchildren sets the number of children to add at each split of the tree
 * @param gcount pointer to a variable, which counts/tracks the number of gates in the tree
 * @param ecount pointer to a variable, which counts/tracks the number of experts in the tree
 * @return Node* pointer to the root node of the newly created tree
 */
Node* Node::createTreeInternal(int depth, int nchildren, int* gcount, int* ecount)
{
    if (depth==0){
        Expert* E=new Expert();
        E->name="E" + std::to_string((*ecount)++);
        return E;
    }
    Gate* root = new Gate();
    root->name="G" + std::to_string((*gcount)++);
    for (int i=0; i<nchildren; i++)
        root->addChild(createTreeInternal(depth-1, nchildren, gcount, ecount));
    return root;
} 



