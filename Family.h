#include <iostream>
#include <vector>

class Family {
    public:
      virtual int is_natural(); // is it a natural link (might not be needed)
      virtual double linkfun(double x); // needs to be implemented by subclass
//    virtual void linkfun_vec(arma::vector<double> x, arma::vector<double> result)
      virtual arma::vector<double> linkfun_vec(arma::vector<double> x); // has a default which applies linkfun to every element, but can be overwritten

}

class ExpertModel {
  // has subclasses NormalModel and GLMModel  (GLMModel uses Family, NormalModel doesn't)
  // has a loglikelhood method
  // has a derivative of loglikelihood method
  // has a sample from conjugate method
   

}
