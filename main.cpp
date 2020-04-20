#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

using namespace std;
using namespace arma;

#include  "NormalExpert.h"

int main(){
vec x;
x<<1<<2<<3;
cout<<x.max()<<endl;
cout<<exp(x).max()<<endl;
return 0;
}

