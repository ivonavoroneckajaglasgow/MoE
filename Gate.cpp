//
// Created by 2022175V on 14/08/2019.
//

#include "Gate.h"
#include "Node.h"
#include "NormalExpert.h"
#include <iostream>
#include <vector>
#include "armadillo"

using namespace arma;
using namespace std;

/**
 * @brief Construct a new Gate:: Gate object.
 * 
 * @param aName a string name of the Gate.
 * @param aParameters a list of Gate parameters for the Gate.
 */
Gate::Gate(string aName,GateParameters aParameters)
        :Node()
{
    this->type="G";
    this->name=aName;
    this->parameters=aParameters;
    cout<<"Gate "<<name<<" has been created."<<endl;
}

Gate::Gate(string aName, GateParameters aParameters,NormalParameters aParameters2, int depth, int nchildren, int* gcount, int* ecount)
        :Node()
{
    this->type="G";
    this->name=aName;
    this->parameters=aParameters;
    this->createTree(aParameters,aParameters2,depth,nchildren,gcount,ecount);
    //this=this->createTree(aParameters,aParameters2,depth,nchildren,gcount,ecount);
    //cout<<"Gate "<<name<<" has been created. It has a depth of "<<depth<<" with "<<nchildren<<" children at each split."<<endl;
} 

Node* Gate::createTree(GateParameters aParameters,NormalParameters aParameters2, int depth, int nchildren, int* gcount, int* ecount)
{
    //Also tried using this instead of root, but then it keeps on overriting in line 46 and can't get the 
    //syntax
    if (depth==0)
        return new NormalExpert("E" + std::to_string((*ecount)++), aParameters2);
    Gate* root = new Gate("G" + std::to_string((*gcount)++), aParameters);
    for (int i=0; i<nchildren; i++)
        root->addChild(createTree(aParameters,aParameters2, depth-1, nchildren, gcount, ecount));
    return root;
} 

/**
 * @brief Destroy the Gate:: Gate object.
 * 
 */
Gate::~Gate(){
    for (int i; i<this->Children.size(); i++)
        delete this->Children[i];
}

/**
 * @brief adds a child to the gate.
 * Adds either an expert or a gate as a child to itself.
 * @param aChild the object to be added as a child.
 */
void Gate::addChild(Node* aChild){
    this -> Children.push_back(aChild);
    aChild -> Parent = this;
    cout<<"Child "<<aChild->name<<" has been added to the parent "<<this->name
        <<" and vice versa."<<endl;
}

/**
 * @brief Prints out the names of all children of the gate.
 * 
 */
void Gate::printChildren(){
    if(Children.size()>1) {
        cout << "Gate " << name << " has " << Children.size() << " children called ";
    }else if(Children.size()==1){
        cout << "Gate " << name << " has a child called ";
    }else{
        cout << "Gate " << name << " does not have children. "<<endl;
    }

    for(int i=0; i<Children.size();i++) {
        cout << Children[i]->name;
        if (i == Children.size() - 2) {
            cout << " and ";
        } else if (i == Children.size()-1) {
            cout << "."<<endl;
        } else {
            cout << ", ";
        }
    }
};

/**
 * @brief Prints out the names of all descendants of the gate.
 * 
 */
void Gate::printDescendants(){

    for(int i=0;i<Children.size();i++){
        cout<<Children[i]->name<<endl;
        if(Children[i]->countChildren()!=0){
            Children[i]->printDescendants();
        }
    }
}

/**
 * @brief Prints out the names of all terminal nodes (experts) descending from the gate.
 * 
 */
void Gate::printTerminalNodes(){

    for (int i = 0; i < Children.size(); i++) {
        if (Children[i]->countChildren() == 0) {
            cout << Children[i]->name << endl;
        } else {
            Children[i]->printTerminalNodes();
        }
    }
}

/**
 * @brief Retrieves all children of the gate.
 * 
 * @return vector<Node*>  a vecor of pointers to all children of the gate.
 */
vector<Node*> Gate::getChildren() {
    return Children;
}

/**
 * @brief An internal function which helps retrieve all descendents of the gate.
 * An internal function that is then called at the Node level to retrieve all descendants of the gate.
 * @param desc vector to be filled in with the descendants of the gate.
 * @return vector<Node*> vector of pointers to the descendants of the gate.
 */
vector<Node*> Gate::getDescendantsInternal(vector<Node*>* desc) {

    for(int i=0; i<this->Children.size();i++){
        desc->push_back(this->Children[i]);
         if(this->Children[i]->countChildren()!=0){
            this->Children[i]->getDescendantsInternal(desc);
        }
    }
    return *desc;
}

/**
 * @brief An internal function which helps retrieve all terminal nodes descending from the gate.
 * An internal function that is then called at the node level to retrieve all terminal nodes of the gate.
 * @param terminal vector to be filled in with the terminal nodes descending from the gate.
 * @return vector<Node*> vector of pointers to the terminal nodes descending from the gate.
 */
vector<Node*> Gate::getTerminalNodesInternal(vector<Node*>* terminal){
    for(int i=0; i<this->Children.size();i++){
        if(this->Children[i]->countChildren()==0){
            cout<<"I can see you are an expert "<<Children[i]->name<<endl;
            terminal->push_back(this->Children[i]);
        }else{
            cout<<"I can see you are a gate "<<Children[i]->name<<endl;
            this->Children[i]->getTerminalNodesInternal(terminal);
        }
    }
    return *terminal;
};

/**
 * @brief An internal function which helps to construct a numeric vector describing a tree.
 * This function is called at the node level and returns a vector containing the number of children 
 * of the gate. Then, the same is performed for each of the children of the gate. If the child is an expert,
 * the describeTreeInternal function will be called from the expert level. 
 * @param description vector to be filled with integers describing a tree.
 * @return vector<int> pointer to a vector of integers with the number of children of the gate.
 */
vector<int> Gate::describeTreeInternal(vector<int>* description){

        description->push_back(Children.size());

        for(int i=0;i<Children.size();i++){
            Children[i]->describeTreeInternal(description);
        }

        return *description;
};

/**
 * @brief Counts the number of children of the gate.
 * 
 * @return int integer number of children of the gate.
 */
int Gate::countChildren(){
    int n;
    n=Children.size();
    return n;
};

/**
 * @brief Counts the number of descendants of the gate.
 * 
 * @return int integer number of descendants of the gate.
 */
int Gate::countDescendants(){
    vector<Node*> desc;
    desc=this->getDescendants();
    return desc.size();
};

/**
 * @brief A top layer function for assigning ID.
 * 
 * Sets the gate id to start at 2 (the root gate is automatically determined in issueID_helper2()) and 
 * the expert id to start at 1. Calls the helper functions to perform the task.
 * 
 */
void Gate::issueID(){

    int gateid=2;
    int expertid=1;

    this->issueID_helper2(&gateid,&expertid);
}
/**
 * @brief First internal function which helps to assign the id's to the nodes of tree vertically.
 * This function checks each child of the Node. If a child doesn't have any children, it assumes that it is an expert and 
 * assigns expert id to it. If a child has some children, a gate id is assigned. 
 * @param gate_id pointer to an integer which tracks gate id's.
 * @param expert_id pointer to an integer which tracks expert id's.
 */
void Gate::issueID_helper1(int* gate_id, int* expert_id){

    for(int i=0;i<this->Children.size();i++){
        if(this->Children[i]->countChildren()==0){
            cout<<"I know "<< this->Children[i]->name <<" is an expert so I assign the ID "<<*expert_id<<endl;
            this->Children[i]->id=*expert_id;
            (*expert_id)++;
        }else{
            cout<<"I know "<< this->Children[i]->name <<" is a gate so I assign ID "<<*gate_id<<endl;
            this->Children[i]->id=*gate_id;
            (*gate_id)++;
        }
    }
}

/**
 * @brief Second internal function which helps to assign the id's to the nodes of tree vertically.
 * This function identifies a root gate by checking the presence of the parent node. If there is no parent 
 * node, an id of 1 is assigned. Next, issueID_helper1() is called to assign the id's to the children of the 
 * gate. The process is repeated for every node in the tree by calling this function recursively.
 * @param gate_id pointer to an integer which tracks gate id's.
 * @param expert_id pointer to an integer which tracks expert id's.
 */
void Gate::issueID_helper2(int* gate_id, int* expert_id){

    if(this->Parent==NULL){
        this->id=1;
        cout<<"I know "<<this->name<<" is a root gate, so I assign ID "<<this->id<<endl;
    }

    if(this->id==1) {
       this->issueID_helper1(gate_id, expert_id);
    }

     for(int i=0;i<this->Children.size();i++){
        this->Children[i]->issueID_helper1(gate_id,expert_id);
     }

    for(int i=0;i<this->Children.size();i++) {
        if (this->Children[i]->countChildren()!=0)
            this->Children[i]->issueID_helper2(gate_id,expert_id);
    }
}


