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
