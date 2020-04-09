//
// Created by 2022175V on 14/08/2019.
//

#include "Node.h"
#include "Gate.h"
#include "NormalExpert.h"

using namespace std;
using namespace arma;


/**
 * @brief Construct a new Node:: Node object.
 * 
 */
Node::Node(){
    this->Parent=NULL;
    this->id=NULL;
}

/**
 * @brief A wrapper for createTreeInternal.
 * Initialises gcount and ecount automatically.  
 * @param aParameters list og parameters for gates.
 * @param aParameters2 list of parameters for experts.
 * @param depth sets the depth of the tree to be created.
 * @param nchildren sets the number of children to add at each split of the tree.
 * @return Node* pointer to the root node of the newly created tree.
 */

Node* Node::createTree(GateParameters aParameters,NormalParameters aParameters2, int depth, int nchildren){
    int gcount=0;
    int ecount=0;
    Node* root;
    root=this->createTreeInternal(aParameters,aParameters2,depth,nchildren,&gcount,&ecount);
    return root;
}
/**
 * @brief Creates a tree object given a set of instructions.
 * Creates a node, which can be either an expert (if depth is zero) or a gate. 
 * @param aParameters list og parameters for gates.
 * @param aParameters2 list of parameters for experts.
 * @param depth sets the depth of the tree to be created.
 * @param nchildren sets the number of children to add at each split of the tree.
 * @param gcount pointer to a variable, which counts/tracks the number of gates in the tree.
 * @param ecount pointer to a variable, which counts/tracks the number of experts in the tree.
 * @return Node* pointer to the root node of the newly created tree.
 */
Node* Node::createTreeInternal(GateParameters aParameters,NormalParameters aParameters2, int depth, int nchildren, int* gcount, int* ecount)
{
    if (depth==0)
        return new NormalExpert("E" + std::to_string((*ecount)++), aParameters2);
    Gate* root = new Gate("G" + std::to_string((*gcount)++), aParameters);
    for (int i=0; i<nchildren; i++)
        root->addChild(createTreeInternal(aParameters,aParameters2, depth-1, nchildren, gcount, ecount));
    return root;
} 

/**
 * @brief Turns a numerical vector describing a tree into a tree.
 * Turns a numerical vector of integers (as produced by Gate->describeTree()) into a tree object.
 * @param description vector of integers containing the number of children of each node.
 * @param aParameters list og parameters for gates.
 * @param aParameters2 list of parameters for experts.
 * @return Node* pointer to the root node of the newly created tree.
 */
Node* Node::translateTree(vector<int> description, GateParameters aParameters, NormalParameters aParameters2) {
    if(accumulate(description.begin(), description.end(), 0)!=description.size()-1)
        cout<<"Warning: the description vector is not complete. All missing entries will be replaced by 1, i.e. an Expert will be added."<<endl;
      // throw "Warning: the description vector is not complete. All missing entries will be replaced by 1, i.e. an Expert will be added.";
    int ecount=0;
    int gcount=0;
    if (description.size()==0)
        return 0;
    if (description[0]==0)
        return new NormalExpert("E" + std::to_string(ecount++), aParameters2);
    Gate* root = new Gate("G" + std::to_string(gcount++), aParameters);
    this->populateGate(root, description, 0, &gcount, &ecount, aParameters, aParameters2);
    return root;
}

/**
 * @brief Internal function, used in the translateTree function.
 * Populates the supplied gate with children. Decides on whether to add a Gate or an Expert as a child based
 * on the supplied description vector.
 * @param parent the gate to be populated.
 * @param description vector of integers containing the number of children of each node.  
 * @param start integer that denotes the position of the description vector.
 * @param gcount pointer to a variable, which counts/tracks the number of gates in the tree.
 * @param ecount pointer to a variable, which counts/tracks the number of experts in the tree. 
 * @param aParameters list og parameters for gates.
 * @param aParameters2 list of parameters for experts.
 * @return int returns the integer denoting last position in the description vector.
 */
int Node::populateGate(Gate* parent, vector<int> description, int start, int* gcount, int* ecount, GateParameters aParameters, NormalParameters aParameters2) {
    int pos =  start;
    for (int i=0; i<description[start]; i++) {
        pos++;
        if (pos >= description.size() || description[pos] == 0) {
            parent->addChild(new NormalExpert("E" + std::to_string((*ecount)++), aParameters2));
        } else {
            Gate *g = new Gate("G" + std::to_string((*gcount)++), aParameters);
            parent->addChild(g);
            pos = populateGate(g, description, pos, gcount, ecount, aParameters, aParameters2);
        }
    }
    return pos;
}

/**
 * @brief Prints the name of the parent of the node.
 * 
 */
void Node::printParent() {
    if (Parent == NULL) {
        if(countChildren()!=0) {
            cout << name << " is a root Gate." << endl;
        }
    } else {
        cout << name << " parent is " << Parent->name << "." << endl;
    }
}

/**
 * @brief Prints names of children of the node.
 * 
 */
void Node::printChildren() {
}

/**
 * @brief Prints named of ancestors of the node.
 * 
 */
void Node::printAncestors() {
    if(Parent==NULL){
        cout <<" " << endl;
    }else{
        cout<<Parent->name<<endl;
        Parent->printAncestors();
    }

}
/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 */
void Node::printDescendants(){

}

/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 */
void Node::printTerminalNodes(){

}

/**
 * @brief Retrieves parent of the node.
 * 
 * @return Gate* pointer to the parent node.
 */
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

}

/**
 * @brief  A dummy function, which is handled at expert and node levels.
 * @return int 
 */
int Node::countChildren() {

}

/**
 * @brief Function which constructs a numeric vector describing a tree.
 * This function returns a vector containing the number of children of each node in the tree. 
 * To do so, the describeTreeInternal() is called on expert and gate level.
 * @param description_test vector to be filled with integers describing a tree.
 * @return vector<int> pointer to a vector of integers with the number of children of the gate.
 */
vector<int> Node::describeTree(){
    vector <int> description_test;
    description_test=this->describeTreeInternal(&description_test);
    return description_test;
}

/**
 * @brief A dummy function, which is handled at expert and node levels.
 * 
 * @param description 
 * @return vector<int> 
 */
vector<int> Node::describeTreeInternal(vector<int>* description){

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
