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
    if (depth==0)
        return new NormalExpert("E" + std::to_string((*ecount)++), aParameters2);
    Gate* root = new Gate("G" + std::to_string((*gcount)++), aParameters);
    for (int i=0; i<nchildren; i++)
        root->addChild(createTree(aParameters,aParameters2, depth-1, nchildren, gcount, ecount));
    return root;
} 

Gate::~Gate(){
    for (int i; i<this->Children.size(); i++)
        delete this->Children[i];
}

void Gate::addChild(Node* aChild){
    this -> Children.push_back(aChild);
    aChild -> Parent = this;
    cout<<"Child "<<aChild->name<<" has been added to the parent "<<this->name
        <<" and vice versa."<<endl;
}

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

void Gate::printDescendants(){

    for(int i=0;i<Children.size();i++){
        cout<<Children[i]->name<<endl;
        if(Children[i]->countChildren()!=0){
            Children[i]->printDescendants();
        }
    }
}

void Gate::printTerminalNodes(){

    for (int i = 0; i < Children.size(); i++) {
        if (Children[i]->countChildren() == 0) {
            cout << Children[i]->name << endl;
        } else {
            Children[i]->printTerminalNodes();
        }
    }
}

vector<Node*> Gate::getChildren() {
    return Children;
}

vector<Node*> Gate::getDescendantsInternal(vector<Node*>* desc) {

    for(int i=0; i<this->Children.size();i++){
        desc->push_back(this->Children[i]);
         if(this->Children[i]->countChildren()!=0){
            this->Children[i]->getDescendantsInternal(desc);
        }
    }
    return *desc;
}

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

vector<int> Gate::describeTreeInternal(vector<int>* description){

        description->push_back(Children.size());

        for(int i=0;i<Children.size();i++){
            Children[i]->describeTreeInternal(description);
        }

        return *description;
};

int Gate::countChildren(){
    int n;
    n=Children.size();
    return n;
};

int Gate::countDescendants(){
    vector<Node*> desc;
    desc=this->getDescendants();
    return desc.size();
};

void Gate::issueID(){

    int gateid=2;
    int expertid=1;

    this->issueID_helper2(&gateid,&expertid);
}

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


