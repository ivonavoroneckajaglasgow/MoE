#include <iostream>
#include <vector>

#include "Family.h"

class PoissonFamily: public Family{
public:
    double linkfun(double mu);
};