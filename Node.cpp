#define _USE_MATH_DEFINES
#define EPS 2.22e-16

#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

#include "Node.h"
#include "Gate.h"

void Node::printParent() {
        cout << name << " parent is " << Parent->name << "." << endl;
}

void Node::printDescendants(){

}

void Node::printTerminalNodes(){

}

Gate* Node::getParent(){
    return Parent;
}

/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 * @return vector<Node*> empty vector.
 */
vector<Node*> Node::getChildren() {
    vector<Node*> empty;
    return empty;
}

/**
 * @brief Retrieves ancestors of the node.
 * Calls the internal function getAncestorsInternal().
 * @return vector<Node*> vector of pointers to the ancestors of the node.
 */
vector<Node*> Node::getAncestors() {
    vector <Node*> ancest_test;
    ancest_test=this->getAncestorsInternal(&ancest_test);
    return ancest_test;
}

/**
 * @brief An internal helper function tht retrieves ancestors of the node.
 * 
 * @param ancest vector to be filled in with ancestors of the node.
 * @return vector<Node*> vector of pointers to the ancestors of the node.
 */
vector<Node*> Node::getAncestorsInternal(vector<Node*>* ancest) {

    if(this->Parent!=NULL){
        ancest->push_back(this->Parent);
        this->Parent->getAncestorsInternal(ancest);
    }
    return *ancest;

}

/**
 * @brief Retrieves the descendents of the node.
 * Calls the getDescendantsInternal() function from gate and expert levels.
 * @return vector<Node*> vector of pointers to the descendants of the node.
 */
vector<Node*> Node::getDescendants() {
    vector <Node*> desc_test;
    desc_test=this->getDescendantsInternal(&desc_test);
    return desc_test;
}

/**
 * @brief A dummy function, which is handled at expert and node levels.
 * @param desc 
 * @return vector<Node*> 
 */
vector<Node*> Node::getDescendantsInternal(vector<Node*>* desc){
vector<Node*> a;
return a;
}

/**
 * @brief Retrieves the terminal nodes descending from the node.
 * Calls the getTerminalNodesInternal() function from gate and expert levels.
 * @return vector<Node*> vector of pointers to the terminal nodes descending from the node.
 */
vector<Node*> Node::getTerminalNodes() {
    vector <Node*> terminal_test;
    terminal_test=this->getTerminalNodesInternal(&terminal_test);
    return terminal_test;
}

/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 * @param terminal 
 * @return vector<Node*> 
 */
vector<Node*> Node::getTerminalNodesInternal(vector<Node*>* terminal){
vector<Node*> a;
return a;
}

int Node::countChildren() {
return 0;
}

/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 */
void Node::issueID() {}
/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 * @param gate_id 
 * @param expert_id 
 */
void Node::issueID_helper1(int* gate_id, int* expert_id){}
/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 * @param gate_id 
 * @param expert_id 
 */
void Node::issueID_helper2(int* gate_id, int* expert_id) {}

int Node::issueIDLR(int start){
    return 0;
}

int Node::leftMostNodeID(){
    return this->idLR;
}

int Node::rightMostNodeID(){
    return 0;
}

int Node::isInRange(Node* node){
    return 0;
}

vec Node::getDescendantRange(){
    vec result;
    result<<this->leftMostNodeID()<<this->rightMostNodeID();
    return result;
}

int Node::countPoints(vector<Node*> node){
   int count=0;
   for(int i=0;i<node.size();i++){
    if(this->isInRange(node[i])==1){
       count=count+1;
    }
   }
   return count;
}

vec Node::getPointIndices(vector<Node*> node){
    vec range=this->getDescendantRange();    
    vec result(this->countPoints(node));
    int position=0;
    for(int i=0;i<node.size();i++){
        if(this->isInRange(node[i])==1){
            result[position]=i;
            position=position+1;
        }
    }
    return result;
}