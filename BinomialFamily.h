#include <iostream>
#include <vector>

#include "Family.h"

class BinomialFamily: public Family{
public:
    double linkfun(double mu);
};