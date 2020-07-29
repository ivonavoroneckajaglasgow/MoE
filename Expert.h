#ifndef MOE_EXPERT_H
#define MOE_EXPERT_H

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#include "ExpertModel.h"
#include "Node.h"

using namespace std;
using namespace arma;

//This is a superclass of NormalExert and ExpertModel objects 
//Most functions are virtual and overwritten at subclass levels

class Expert: public Node{
    public:
    vec y;       //data observed
    vec eta;     // XB
    mat X;
    vec beta;
    double var;  // some variance parameter  
    ExpertModel* expertmodel; // the model for this expert
    Expert();//constructor
};

#endif //MOE_EXPERT_H