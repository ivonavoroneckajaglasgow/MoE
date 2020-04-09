#include <iostream>
#include <vector>

#include "Family.h"

class ExpertModel {
    public:
    Family  family;
    double  loglik(double x); //NEEDED: some sort of way of feeding in parameters generically for all models
    double  dloglik(double x);
    double  sampleConjugate(double x);
};