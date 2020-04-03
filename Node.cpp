//
// Created by 2022175V on 14/08/2019.
//

#include "Node.h"
#include "Gate.h"
#include "NormalExpert.h"

using namespace std;

Node::Node(){
    this->Parent=NULL;
    this->id=NULL;
}

void Node::printParent() {
    if (Parent == NULL) {
        if(countChildren()!=0) {
            cout << name << " is a root Gate." << endl;
        }
    } else {
        cout << name << " parent is " << Parent->name << "." << endl;
    }
}

void Node::printChildren() {
}

void Node::printAncestors() {
    if(Parent==NULL){
        cout <<" " << endl;
    }else{
        cout<<Parent->name<<endl;
        Parent->printAncestors();
    }

}

void Node::printDescendants(){

};

void Node::printTerminalNodes(){

}

Gate* Node::getParent(){
    return Parent;
}

vector<Node*> Node::getChildren() {
    vector<Node*> empty;
    return empty;
}

vector<Node*> Node::getAncestors() {
    vector <Node*> ancest_test;
    ancest_test=this->getAncestorsInternal(&ancest_test);
    return ancest_test;
}

vector<Node*> Node::getAncestorsInternal(vector<Node*>* ancest) {

    if(this->Parent!=NULL){
        ancest->push_back(this->Parent);
        this->Parent->getAncestorsInternal(ancest);
    }
    return *ancest;

}

vector<Node*> Node::getDescendants() {
    vector <Node*> desc_test;
    desc_test=this->getDescendantsInternal(&desc_test);
    return desc_test;
}

vector<Node*> Node::getDescendantsInternal(vector<Node*>* desc){

};

vector<Node*> Node::getTerminalNodes() {
    vector <Node*> terminal_test;
    terminal_test=this->getTerminalNodesInternal(&terminal_test);
    return terminal_test;
}

vector<Node*> Node::getTerminalNodesInternal(vector<Node*>* terminal){

};

int Node::countChildren() {

}

vector<int> Node::describeTree(){
    vector <int> description_test;
    description_test=this->describeTreeInternal(&description_test);
    return description_test;
}

vector<int> Node::describeTreeInternal(vector<int>* description){

}

void Node::issueID() {}

void Node::issueID_helper1(int* gate_id, int* expert_id){}

void Node::issueID_helper2(int* gate_id, int* expert_id) {}
