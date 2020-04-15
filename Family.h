#include <iostream>
#include <vector>


using namespace std;
//using namespace arma;

class Family {
    public:
      //virtual int is_natural(); // is it a natural link (might not be needed)
      virtual double linkfun(double mu); // needs to be implemented by subclass
      vector<double> linkfun_vec(vector<double> x);
      //virtual arma::vector<double> linkfun_vec(arma::vector<double> x); // has a default which applies linkfun to every element, but can be overwritten
};
