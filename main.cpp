#define _USE_MATH_DEFINES
#define EPS 1e-16

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"
#include <chrono>

using namespace std;
using namespace arma;
using namespace std::chrono; 

#include "Gate.h"

int main(){
cout<<"SETTING UP GATE"<<endl;

mat X={{0.1,0.4,0.7,1},
       {0.2,0.5,0.8,1.1},
       {0.3,0.6,0.9,1.2}};
X.print("Design matrix X:");
int n=X.n_rows;
cout<<"Number of observations n="<<n<<endl;
int p=X.n_cols;
cout<<"Number of covariates p="<<p<<endl;

mat z={{0,1},
       {0,0},
       {1,0}};

int r=z.n_cols;
cout<<"Number of splits (excluding ref class) r="<<r<<endl;
z.print("Matrix of allocations z:");
cout<<"Each row corresponds to an observation and each column to a split."<<endl;

Gate* G= new Gate();

vec gamma("-0.84085548,1.38435934,-1.25549186,0.07014277,1.71144087,-0.60290798,-0.47216639,-0.63537131");
gamma.print("Gating parameters gamma:");
cout<<"Length of gamma (rp): "<<gamma.size()<<endl;

mat pi=G->pi_calculator(X,gamma);
pi.print("pi calculator 1:");

mat pi2=G->pi_calculator2(X,gamma);
pi2.print("pi calculator 2:");

mat pi3=G->pi_calculator2(X,gamma);
pi3.print("pi calculator 3:");

cout<<"Simulation:"<<endl;

int N=10;

mat runtime(N,2);

int n_star=10000;
int p_star=10000;
int r_star=5;

cout<<"n="<<n_star<<endl;
cout<<"p="<<p_star<<endl;
cout<<"r="<<r_star<<endl;

for(int i=0; i<N; i++){
cout<<i<<endl;
mat X_star(n_star,p_star,fill::randu);
vec gamma_star(r_star*p_star,fill::randn);

auto start= high_resolution_clock::now();
clock_t c_start = clock();
mat pi=G->pi_calculator(X_star,gamma_star);
//mat pi=G->pi_calculator2(X_star,gamma_star);
//mat pi=G->pi_calculator3(X_star,gamma_star);
clock_t c_end = clock();
auto stop = high_resolution_clock::now(); 
auto duration = duration_cast<microseconds>(stop - start); 
double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;

runtime(i,0)=time_elapsed_ms;
runtime(i,1)=duration.count()/1000;
}

cout<<"Pi calculator:"<<endl;
cout<<"Average CPU and System time used (respectively):"<<endl;
cout<<mean(runtime,0)<<" ms\n";

return 0; 
}