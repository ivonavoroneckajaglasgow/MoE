#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <cmath>
#include "armadillo"

#define _USE_MATH_DEFINES
#define EPS 2.22e-16
#define SigmaMultiple 2

#include "Gate.h"
#include "Expert.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Gate:: Gate object
 * 
 */
Gate::Gate(){
    //cout<<"Gate has been created."<<endl;
}

/**
 * @brief Destroy the Gate:: Gate object
 * 
 */
Gate::~Gate(){
    cout<<"Gate has been deleted."<<endl;
}


/**
 * @brief Construct a new Gate:: Gate object
 * Copy constructor
 * @param gate 
 */
Gate::Gate(const Gate &gate){
    cout<<"A copy of Gate has been created."<<endl;
}

/**
 * @brief Creates a new gate that is a copy of this one
 * 
 * @return Gate::Gate* pointer to the new gate 
 */
 Gate* Gate::copyThis(){
     //cout<<"Making a copy of "<< name<<endl;
     Gate* Copy=new Gate();
     Copy->name=this->name;
     Copy->gamma=this->gamma;
     Copy->Omega=this->Omega;
     Copy->idLR=this->idLR;
     return Copy;
 }
 
/**
 * @brief Creates a copy of the gate and all its descedants and stores that in Copy
 * 
 * @param Copy pointer to the copy of the gate with its descendants
 */
void Gate::copyStructure(Gate* Copy){ 
    for(int i=0; i<this->countChildren();i++){
         if(this->Children[i]->countChildren()==0){
             Copy->addChild(this->Children[i]->copyThis());
         }else{
             Copy->addChild(this->Children[i]->copyThis());
             //cout<<"I see that "<<Copy->Children[i]->name<<" is a gate, so I continue"<<endl;
             dynamic_cast<Gate*>(this->Children[i])->copyStructure(dynamic_cast<Gate*>(Copy->Children[i]));
         }
     }
 }

/**
 * @brief Wrapper for copyStructure(Gate* gate)
 * 
 * @return Gate* pointer to the copied gate with all its descendants
 */
 Gate* Gate::copyStructure(){
     Gate* Copy=this->copyThis();
     this->copyStructure(Copy);
     return Copy;
 }

void Gate::replaceGate(Gate* replacement){
     this->name=replacement->name;
     this->gamma=replacement->gamma;
     this->Omega=replacement->Omega;
     this->Children=replacement->Children;
}

void Gate::replaceChild(int which, Node* newChild){
    cout<<"Replace the "<<which<<"-th child of "<<this->name<<endl;
    cout<<"Replace "<< this->Children[which]->name <<" by "<<newChild->name<<endl;
    cout<<"Set "<<this->name<<" to be the parent of "<<newChild->name<<endl;
    this->Children[which]=newChild;
    newChild->Parent=this;
}

// void Gate::exchangeWith(int which, Node* other){
//     Gate* myParent=this->Parent;
//     Gate* otherParent=other->Parent;
//     myParent->replaceChild(which,other);
//     otherParent->replaceChild(which,this);
// }

/**
 * @brief Add a child to the gate and make gate the parent of child
 * 
 * @param aChild pointer to the node to be added as a child 
 */
void Gate::addChild(Node* aChild){
    this -> Children.push_back(aChild);
    aChild -> Parent = this;
    //cout<<"Child "<<aChild->name<< " has been added to the parent "<<name<<"."<<endl;
}

/**
 * @brief Deletes children and itself as a parent
 * 
 */
void Gate::deleteChildren(){
    for(int i=0;i<this->countChildren();i++){
        this->Children[i]->Parent=NULL;
    }
    vector<Node*> empty;
    this->Children=empty;
}

/**
 * @brief Prints out the names of all children
 * 
 */
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
}

/**
 * @brief Prints out the names of all descendants
 * 
 */
void Gate::printDescendants(){

    for(int i=0;i<Children.size();i++){
        cout<<Children[i]->name<<endl;
        if(Children[i]->Children.size()!=0){
            Children[i]->printDescendants();
        }
    }
}

/**
 * @brief Prints out the names of all terminal nodes (experts)
 * 
 */
void Gate::printTerminalNodes(){
    
    for (int i = 0; i < Children.size(); i++) {
        if (Children[i]->Children.size() == 0) {
            cout << Children[i]->name << endl;
        } else {
            Children[i]->printTerminalNodes();
        }
    }
}

/**
 * @brief Returns all children 
 * 
 * @return vector<Node*> vector of pointers to the children
 */
vector<Node*> Gate::getChildren() {
    return Children;
}

/**
 * @brief An internal function which outputs all descendents 
 * This function is called at the Node level 
 * @param desc vector to be filled in with the descendants 
 * @return vector<Node*> vector of pointers to the descendants 
 */
vector<Node*> Gate::getDescendantsInternal(vector<Node*>* desc) {

    for(int i=0; i<this->Children.size();i++){
        desc->push_back(this->Children[i]);
         if(this->Children[i]->countChildren()!=0){
            this->Children[i]->getDescendantsInternal(desc);
        }
    }
    return *desc;
}

/**
 * @brief An internal function which outputs all terminal nodes 
 * This function is called at the Node level 
 * @param terminal vector to be filled in with the terminal nodes 
 * @return vector<Node*> vector of pointers to the terminal nodes 
 */
vector<Node*> Gate::getTerminalNodesInternal(vector<Node*>* terminal){
    for(int i=0; i<this->Children.size();i++){
        if(this->Children[i]->countChildren()==0){
            //cout<<"I can see you are an expert "<<Children[i]->name<<endl;
            terminal->push_back(this->Children[i]);
        }else{
            //cout<<"I can see you are a gate "<<Children[i]->name<<endl;
            this->Children[i]->getTerminalNodesInternal(terminal);
        }
    }
    return *terminal;
}

/**
 * @brief Returns the total number of children
 * 
 * @return int integer number of children
 */
int Gate::countChildren(){
    return static_cast<int>(Children.size());
}

/**
 * @brief returns the total number of descendants 
 * 
 * @return int integer number of descendants 
 */
int Gate::countDescendants(){
    vector<Node*> desc;
    desc=this->getDescendants();
    return static_cast<int>(desc.size());
}

/**
 * @brief A top layer function for assigning ID
 * 
 * Sets the gate id to start at 2 (the root gate is automatically determined in issueID_helper2() and 
 * the expert id to start at 1. Calls the helper functions to perform the task.
 * 
 */
void Gate::issueID(){

    int gateid=2;
    int expertid=1;

    this->issueID_helper2(&gateid,&expertid);
}

/**
 * @brief First internal function which helps to assign the IDs to the nodes of tree vertically
 * This function checks each child of every node. If a child doesn't have any children, it assumes that it is an expert and 
 * assigns expert id to it. If a child has some children, a gate id is assigned 
 * @param gate_id pointer to an integer which tracks gate IDs
 * @param expert_id pointer to an integer which tracks expert IDs
 */
void Gate::issueID_helper1(int* gate_id, int* expert_id){

    for(int i=0;i<this->Children.size();i++){
        if(this->Children[i]->countChildren()==0){
            //cout<<"I know "<< this->Children[i]->name <<" is an expert so I assign the ID "<<*expert_id<<endl;
            this->Children[i]->id=*expert_id;
            (*expert_id)++;
        }else{
            //cout<<"I know "<< this->Children[i]->name <<" is a gate so I assign ID "<<*gate_id<<endl;
            this->Children[i]->id=*gate_id;
            (*gate_id)++;
        }
    }
}

/**
 * @brief Second internal function which helps to assign the IDs to the nodes of tree vertically
 * This function identifies a root gate by checking the presence of the parent node. If there is no parent 
 * node, an id of 1 is assigned. Next, issueID_helper1() is called to assign the IDs to the children of the 
 * gate. The process is repeated for every node in the tree by calling this function recursively
 * @param gate_id pointer to an integer which tracks gate IDs
 * @param expert_id pointer to an integer which tracks expert IDs
 */
void Gate::issueID_helper2(int* gate_id, int* expert_id){

    if(this->Parent->Parent==NULL){
        this->id=1;
        //cout<<"I know "<<this->name<<" is a root gate, so I assign ID "<<this->id<<endl;
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

/**
 * @brief Log-likelihood function for the model 
 * z~Multinom(pi_1,...,pi_r)
 * @param z matrix of allocations
 * @param pi matrix of mixing proportions
 * @return double log-likelihood value
 */
double Gate::loglik(mat z, mat pi){
    return static_cast<double>(sum(vectorise(z%log(pi)))+sum((1-sum(z,1))%log(1-sum(pi,1))));
}

/**
 * @brief Log-likelihood function for the model 
 * z~Multinom(pi_1,...,pi_r) 
 * @param X design matrix
 * @param gamma vector of gating parameters
 * @param z matrix of allocations
 * @return double log-likelihood value
 */
double Gate::loglik(mat X, vec gamma, mat z){
    int p=static_cast<int>(X.n_cols);
    int r=static_cast<int>(z.n_cols);
    mat gamma2=gamma;
    gamma2.reshape(p,r);
    return sum(sum(z%(X*gamma2),1)-log(1+sum(exp(X*gamma2),1)));
}

/**
 * @brief mixing proportions calculator
 * 
 * @param X design matrix
 * @param gamma vector of gating parameters
 * @return mat matrix of mixing proportions (rows - observations, columns - splits)
 */
mat Gate::pi_calculator(mat X, vec gamma){
    int p=static_cast<int>(X.n_cols);
    int r=static_cast<int>(gamma.size()/p);
    mat gamma2=gamma;
    gamma2.reshape(p,r);
    mat helper=exp(X*gamma2);
    mat rowsums=this->getRowSumsMat(helper);
    mat final=helper/(1+rowsums);
    mat E(final.n_rows,final.n_cols); 
    E.fill(sqrt(EPS));
    mat zeromat(final.n_rows,final.n_cols);
    zeromat.zeros();
    mat pi=arma::max(final,E);
    mat error=arma::max(zeromat,this->getRowSumsMat(pi)+E-1);
    mat LHS = pi-(pi+E)/this->getRowSumsMat(pi+E)%error;
    return arma::max(LHS,E);
}

/**
 * @brief returns rowsums of a matrix A
 * 
 * @param A matrix 
 * @return mat matrix that contains copies of row sums in every column 
 */
mat Gate::getRowSumsMat(mat A){
    mat rowsums(A.n_rows,A.n_cols);
    for(int i=0; i<A.n_cols; i++){
        rowsums.col(i)=sum(A,1);
    }
    return rowsums;
}

/**
 * @brief mixing proportions calculator (loops over each observation)
 * 
 * @param X design matrix
 * @param gamma vector of gating parameters
 * @return mat matrix of mixing proportions (rows - observations, columns - splits)
 */
 mat Gate::pi_calculator2(mat X, vec gamma) {
    int p = static_cast<int>(X.n_cols);
    int rp= static_cast<int>(gamma.size());
    int r = rp/p;
    double current;
    mat result(X.n_rows,r);
    for (int i=0; i<X.n_rows; i++) {
          double sum=1;
          for (int k=0; k<r; k++) {
            current = as_scalar(exp(X.row(i)*gamma(span(k*p, (k+1)*p-1)))); //might want gamma as a matrix
            sum += current;
            result(i,k) = current;
        } 
          for (int k=0; k<r; k++) {
            current = result(i,k) / sum;
            if (current < sqrt(EPS))
            current=sqrt(EPS);
            if (current > 1-sqrt(EPS))
            current=1-sqrt(EPS);
            result(i,k) = current;
         }
    }
    return result;
}       

/**
 * @brief mixing proportions calculator (subtracting max)
 * 
 * @param X design matrix
 * @param gamma vector of gating parameters
 * @return mat matrix of mixing proportions (rows - observations, columns - splits)
 */
mat Gate::pi_calculator3(mat X, vec gamma) {
    int p = static_cast<int>(X.n_cols);
    int rp= static_cast<int>(gamma.size());
    int r = static_cast<int>(rp/p);
    mat result(X.n_rows,r);
    double current;
    // try reshaping gamma here
    for (int i=0; i<X.n_rows; i++) {
        double max=0;
        for (int k=0; k<r; k++) {
            current = as_scalar(exp(X.row(i)*gamma(span(k*p, (k+1)*p-1))));
            if (current>max)
                max=current;
            result(i,k) = current;
        }
        double sum = exp(-max);
        for (int k=0; k<r; k++) {
            result(i,k) = exp(result(i,k)-max);
            sum += result(i,k);
        }
        for (int k=0; k<r; k++) {
            current = result(i,k) / sum;
            if (current < sqrt(EPS))
                current=sqrt(EPS);
            if (current > 1-sqrt(EPS))
                current=1-sqrt(EPS);
            result(i,k) = current;
        }
    }
    return result;
}       

/**
 * @brief Score function used in IWLS
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param pi matrix of mixing proportions
 * @return vec vector of scores
 */
vec Gate::score(mat X, mat z, mat pi){
    return vectorise(X.t()*(z-pi));
}

/**
 * @brief Hessian matrix calculator used in IWLS
 * 
 * @param X design matrix
 * @param pi matrix of mixing proportions
 * @return mat matrix containing Hessian
 */
mat Gate::hessian(mat X, mat pi){
 int r=static_cast<int>(pi.n_cols);
 int p=static_cast<int>(X.n_cols);

 mat H(r*p,r*p);

 for(int i=0; i<r; i++){
     for(int j=0; j<r; j++){
         mat a;
         if(i==j) a=-pi.col(i)%(1-pi.col(i));
         if(i!=j) a= pi.col(i)%pi.col(j);
        mat a2(X.n_rows,X.n_cols);
            for(int i=0; i<X.n_cols; i++){
            a2.col(i)=a;
        }
    H.submat(span(i*p,(i+1)*p-1),span(j*p,(j+1)*p-1))=X.t()*(a2%X); 
   }
}
 
return H;
}

/**
 * @brief Finds the MLE of gating parameters gamma using the Cholesky decomposition
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix of gamma
 * @param L pointer to the L matrix in the Cholesky decomposition
 * @return vec MLE of gamma
 */
vec Gate::findGammaMLEChol(mat X, mat z, mat Omega, mat* L){
    vec gamma(X.n_cols*z.n_cols);
    gamma.zeros();
     for(int i=0; i<100; i++){
        mat pi=this->pi_calculator(X,gamma);
        vec gamma_old=gamma;
        mat LHS=this->hessian(X,pi)-Omega;
        *L=chol(-LHS);
        mat RHS=this->score(X,z,pi)-Omega*gamma;
        vec a=solve((*L).t(),RHS);
        gamma=gamma+solve(*L,a);
    if(all(abs(gamma-gamma_old)<(sqrt(EPS),sqrt(EPS)*abs(gamma)).max())) break;
    }
    return gamma;
}

/**
 * @brief Finds the MLE of gating parameters gamma using the Cholesky decomposition
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix of gamma
 * @return vec MLE of gamma
 */
vec Gate::findGammaMLEChol(mat X, mat z, mat Omega){
    mat L;
    return this->findGammaMLEChol(X,z,Omega,&L);
}

/**
 * @brief Finds the MLE of gating parameters gamma using the QR decomposition
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix of gamma
 * @param R pointer to the R matrix in the QR decomposition
 * @return vec MLE of gamma
 */
vec Gate::findGammaMLEQR(mat X, mat z, mat Omega, mat* R){
    vec gamma(X.n_cols*z.n_cols);
    gamma.zeros();
    for(int i=0; i<100; i++){
        mat pi=this->pi_calculator(X,gamma);
        vec gamma_old=gamma;
        mat LHS=this->hessian(X,pi)-Omega;
        mat RHS=this->score(X,z,pi)-Omega*gamma;
        mat Q;
        qr(Q,*R,LHS);
        gamma=gamma-solve(*R,Q.t()*RHS);
    if(all(abs(gamma-gamma_old)<(sqrt(EPS),sqrt(EPS)*abs(gamma)).max())) break;
    }
    return gamma;
}

/**
 * @brief Finds the MLE of gating parameters gamma using the QR decomposition
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix of gamma
 * @return vec MLE of gamma
 */
vec Gate::findGammaMLEQR(mat X, mat z, mat Omega){
    mat R;
    return this->findGammaMLEQR(X,z,Omega,&R);
}

/**
 * @brief Finds the MLE of gating parameters gamma 
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix of gamma
 * @return vec MLE of gamma
 */
vec Gate::findGammaMLE(mat X, mat z, mat Omega){
    int p=static_cast<int>(X.n_cols);
    int r=static_cast<int>(z.n_cols);
    vec gamma(p*r);
    gamma.zeros();
    vec diagonals(p*r);
    diagonals.fill(0.001);
    Omega=diagmat(diagonals);
    for(int i=0; i<100; i++){
        mat pi=this->pi_calculator(X,gamma);
        vec gamma_old=gamma;
        gamma=gamma-solve(this->hessian(X,pi)-Omega,this->score(X,z,pi)-Omega*gamma);
    if(all(abs(gamma-gamma_old)<(sqrt(EPS),sqrt(EPS)*abs(gamma)).max())) break;
    }
    return gamma;
}

/**
 * @brief Funtion that takes in a vector of mixing proportions for one point and returns the Cholesky decomposion of matrx A used in the IWLS estimation of gamma
 * 
 * @param pi vector of mixing proportions for one point
 * @return mat the Cholesky decomposion of matrx A used in the IWLS estimation of gamma
 */
mat Gate::makeAchol(vec pi){
  mat helper(pi.size(),pi.size());
  for(int i=0; i<pi.size();i++){
      for(int j=0; j<pi.size(); j++){
          helper(i,j)=pi[i]*pi[j];
      }
  }
  return chol(diagmat(pi)-helper);
}

/**
 * @brief Creates new X matrix used in the IWLS estimation of gamma
 * 
 * @param X design matrix
 * @param pi matrix of mixing proportions
 * @return mat new X matrix used in the IWLS estimation of gamma
 */
mat Gate::getXout(mat X, mat pi){
    int r=static_cast<int>(pi.n_cols);
    int n=static_cast<int>(X.n_rows);
    int p=static_cast<int>(X.n_cols);
    mat Xout(n*r,p*r);

    for(int j=0; j<n; j++){
        mat C=this->makeAchol(vectorise(pi.row(j)));
        mat helper(p*C.n_cols,C.n_cols);
        for(int i=0; i<p; i++){
             vec helper1=vectorise(X.row(j));
             helper.rows(i*r,(i+1)*r-1)=helper1[i]*C;
        }
        helper.reshape(helper.n_cols,helper.n_rows);
        Xout.rows(j*r,(j+1)*r-1)=helper;
    }

    return Xout;
}

/**
 * @brief Creates the new response used in the IWLS estimation of gamma
 * 
 * @param z matrix of allocations
 * @param pi matrix of mixing proportions
 * @return vec new response
 */
vec Gate::getZeta(mat z, mat pi){
    int r=static_cast<int>(pi.n_cols);
    int n=static_cast<int>(pi.n_rows);
    vec zeta(n*r);
    for(int i=0; i<n; i++){
        mat C=this->makeAchol(vectorise(pi.row(i)));
        mat helper=z.row(i)-pi.row(i);
        helper.reshape(helper.n_cols,helper.n_rows);
        zeta(span(i*r,(i+1)*r-1))=solve(C.t(),helper);
}
return zeta;
}

/**
 * @brief Finds the estimate of gamma using IWLS and the QR decomposition
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix for gamma
 * @param R pointer to matrix R in the QR decomposition
 * @return vec estimate of gamma
 */
vec Gate::findGammaQR(mat X, mat z, mat Omega, mat* R){
    int r=static_cast<int>(z.n_cols);
    int n=static_cast<int>(X.n_rows);
    int p=static_cast<int>(X.n_cols);
    vec diagonals(p*r);
    diagonals.fill(0.001);
    Omega=diagmat(diagonals);
    mat Omega_chol=chol(Omega);
    vec gamma(p*r);
    gamma.zeros();
    if(n==0){
        *R=Omega_chol;
        return gamma;
    }
    for(int i=0; i<100; i++){
        vec gamma_old=gamma;
        mat pi=this->pi_calculator(X,gamma);
        mat Xout=this->getXout(X,pi);
        vec zeta=this->getZeta(z,pi);
        mat LHS=join_cols(Xout,Omega_chol);
        mat RHS(LHS.n_rows,1);
        RHS.zeros();
        RHS.rows(0,zeta.size()-1)=Xout*gamma+zeta;
        mat Q;
        qr(Q,*R,LHS);
        gamma=solve(*R,Q.t()*RHS);
        if(all(abs(gamma-gamma_old)<(sqrt(EPS),sqrt(EPS)*abs(gamma)).max())) break;
    }
    return gamma;
}

/**
 * @brief Finds the estimate of gamma using IWLS and the QR decomposition
 * 
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix for gamma
 * @return vec estimate of gamma
 */
vec Gate::findGammaQR(mat X, mat z, mat Omega){
    mat R;
    return this->findGammaQR(X,z,Omega,&R);
}

/**
 * @brief Function proposes a new value of gamma and checks its acceptance. If accepted returns the new value of gamma, if rejected the old value
 * 
 * @param gammaold current gamma
 * @param X design matrix
 * @param z matrix of allocations
 * @param Omega prior variance-covariance matrix of gamma
 * @return vec new gamma
 */
vec Gate::updateGamma(vec gammaold, mat X, mat z, mat Omega){
   vec diagonals(gammaold.size());
   diagonals.fill(0.001);
   Omega=diagmat(diagonals);
   vec mu_gamma(gammaold.size());
   mu_gamma.zeros();
   mat R;
   vec gammahat = this->findGammaQR(X, z, Omega,&R);
   vec v(gammaold.size(),fill::randn);
   mat Sigma_gamma=Omega.i(); 
   mat RHS(R.n_rows,1);
   RHS.zeros();
   RHS.rows(0,gammahat.size()-1)=v;
   vec gammanew=gammahat+sqrt(SigmaMultiple)*solve(R,RHS); //random mvn
   double loglik_old=this->loglik(z,this->pi_calculator(X,gammaold));
   double loglik_new=this->loglik(z,this->pi_calculator(X,gammanew));
   double proposal_old=sum(this->logmvndensity(gammaold,gammahat,&R));
   double proposal_new=sum(this->logmvndensity(gammanew,gammahat,&R));
   double prior_old=sum(this->logmvndensity(gammaold,mu_gamma, Sigma_gamma));
   double prior_new=sum(this->logmvndensity(gammanew,mu_gamma, Sigma_gamma));
   double acceptance=loglik_new-loglik_old+proposal_old-proposal_new+prior_new-prior_old;
   double u=randu();
   bool accept=u<exp(acceptance);
  if(accept==1){
      return gammanew;
  }else{
  return gammaold; 
  }
}

/**
 * @brief Mutivariate normal density on a log scale
 * 
 * @param response response vector
 * @param mean mean vector
 * @param Sigma variance covariance matrix
 * @return vec mutivariate normal density on a log scale
 * 
 */
vec Gate::logmvndensity(vec response, vec mean, mat Sigma){
   int k = static_cast<int>(Sigma.n_cols);
   //return 1/(pow(2*M_PI,k/2)*sqrt(det(Sigma)))*exp(-0.5*(response-mean).t()*Sigma.i()*(response-mean)); - not log scale
   return -k/2*log(2*M_PI)-0.5*log(det(Sigma))-0.5*(response-mean).t()*Sigma.i()*(response-mean);
}

/**
 * @brief Mutivariate normal density on a log scale
 * 
 * @param response response vector
 * @param mean mean vector
 * @param R pointer to the R matrix from QR decomposition
 * @return vec mutivariate normal density on a log scale
 * 
 */
vec Gate::logmvndensity(vec response, vec mean, mat* R){
int k=static_cast<int>((*R).n_rows);
//return -k/2*log(2*M_PI)+0.5*sum(log(pow((*R).diag(),2)))-0.5*(response-mean).t()*((*R).t()*(*R))*(response-mean);
mat result;
result =-k/2*log(2*M_PI)+sum(log((*R).diag()))-0.5*sum(pow((*R)*(response-mean),2));
return vectorise(result);
}

/**
 * @brief Issues IDs to nodes left to right
 * 
 * @param start starting index (always zero and wrapped in the net function)
 * @return int finishing index
 */
int Gate::issueIDLR(int start){
    this->idLR=start++;
    for(int i=0;i<countChildren();i++){
        start=this->Children[i]->issueIDLR(start);
    }
    return start;
}

/**
 * @brief Issues IDs to nodes left to right
 * 
 * @return int recursive function
 */
int Gate::issueIDLR(){
    int start=0;
    return this->issueIDLR(start);
}

/**
 * @brief Finds the ID of the righ most node
 * 
 * @return int the ID of the right most node
 */
 int Gate::rightMostNodeID(){
     return this->Children[countChildren()-1]->rightMostNodeID();
 }

/**
 * @brief Checks if node is in the range of descendants
 * 
 * @param node node to check
 * @return int yes or no 
 */
int Gate::isInRange(Node* node){
    vec range=this->getDescendantRange(); 
    return range[0] <= node->idLR && node->idLR <= range[1];
}

/**
 * @brief Returns which child of gate leads to node
 * 
 * @param node node of interest
 * @return int child index
 */
int Gate::whichSide(Node* node){
    int check;
    for(int i=0; i<this->countChildren();i++){
        if(this->Children[i]->countChildren()!=0){
            check=this->Children[i]->isInRange(node);
            if(check==1){
                return i;
            }
        }
        if(this->Children[i]->idLR==node->idLR){
            cout<<"Node is my child"<<endl;
            return this->whichChild(node);
        }
    }
}

/**
 * @brief returns a vector of length of number of splits, entry of 1 indicates which split node is in
 * 
 * @param node node to check
 * @return rowvec vector where 1 indicates which split node is in
 */
rowvec Gate::getZ_perpoint(Node* node){
    rowvec z(this->countChildren());   
    int    test=this->isInRange(node); //check if node is in the list of descendants for this gate
    if(test==1){ //if it is
        for(int i=0;i<this->countChildren();i++){ //for each child/split of the gate
            z[i]=this->Children[i]->isInRange(node);//check which split this point went down
        }
    }else{
        z.fill(-1); // if node is not in the list of descendants, then return -1
    }
return z;
}

/**
 * @brief returns allocations matrix z
 * 
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 * @return mat allocations matrix z
 */
mat Gate::getZ(vector<Node*> z_final){
mat z(1,this->countChildren());
z.fill(-1);
for(int i=0;i<z_final.size();i++){
rowvec v=this->getZ_perpoint(z_final[i]);
if(sum(v)>=0) z=join_cols(z,v); //if that point is relevant (not -1) add it to the final z matrix
}
z.shed_row(0); //get rid of the first row
z.shed_col(0); //get rid of the reference column
return z;
}

/**
 * @brief updates allocations z
 * 
 * @param y response vector
 * @param X design matrix
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 * @return vector<Node*> updated allocations vector of pointers
 */
vector<Node*> Gate::updateZ(vec y, mat X,vector<Node*> z_final){
    vector<Node*> result=z_final;
    //IF ROOT UPDATE ALL TO BE ADDED AND FIXED
    // if(this->Parent==NULL){
    //     for(int i=0;i<X.n_rows;i++){
    //     vec alpha=this->getSampleProbsForZ(y,X.row(i));
    //     result[i]=this->updateZ_onepoint_sample(alpha);
    // }
    // }else{
    vec points=this->getPointIndices(z_final);
    for(int i=0;i<points.size();i++){
        vec index(1);
        index.fill(points[i]);
        int current=static_cast<int>(points[i]);
        vec alpha=this->getSampleProbsForZ(this->subsetY(y,index),X.row(current));
        result[current]=this->updateZ_onepoint_sample(alpha);
    }
    //} 
    return result;
}

/**
 * @brief Produces the vector of probabilities with which one point should be assigned to terminal nodes (experts)
 * 
 * @param y response (one point)
 * @param X relevant row of the design matrix (one point)
 * @return vec vector of probabilities
 */
vec Gate::getSampleProbsForZ(vec y, mat X){ //y and X are for one i
     vector<Node*> terminals=this->getTerminalNodes();
     vec alpha(terminals.size());
     for(int i=0;i<terminals.size();i++){
         Expert* current=dynamic_cast<Expert*>(terminals[i]);
         vec beta=current->beta;
         double logsigma_sq=current->logsigma_sq;
         alpha[i]=log(this->getPathProb(current,X))+current->expertmodel->loglik(y,X*beta,logsigma_sq);
        }
     double sums=sum(alpha);
     alpha=alpha/sums;
     return alpha;
 }

/**
 * @brief Given a vector of probabilities, assigns one point to one of the terminal nodes 
 * 
 * @param alpha vector of probabilities
 * @return Node* pointer to the node which point is assigned to
 */
 Node* Gate::updateZ_onepoint_sample(vec alpha){
     if(sum(alpha)!=1){
     double sums=sum(alpha);
     alpha=alpha/sums;   
     }
     vector<Node*> terminals=this->getTerminalNodes();
     double rnum = randu();
        for(int i=0; i<alpha.size(); i++){
            if (rnum < alpha[i]) return terminals[i];
            else{
                rnum -= alpha[i];
            }
        }
        cout<<"we should never be here update Z"<<endl;
        alpha.print("alphas:");
        cout<<"rnum "<<rnum<<endl;
    return terminals[alpha.size()-1];
 }

/**
 * @brief Returns a vector of children's left to right IDs 
 * 
 * @return vec vector of children's left to right IDs 
 */
 vec Gate::getChildrenIndicesLR(){
     vec result(this->countChildren());
     for(int i=0; i<this->countChildren(); i++)
     result[i]=this->Children[i]->idLR;
     return result;
 }

/**
 * @brief returns an integer indicating which child is node
 * used to check whether node is in the reference split or not
 * @param node node to check
 * @return int integer indicating which child is node, if not a child returns -1
 */
 int Gate::whichChild(Node* node){
     vec children=this->getChildrenIndicesLR();
     uvec index=find(children==node->idLR);
     if(index.n_rows==0){
        return -1;
     }else{
     return static_cast<int>(as_scalar(index));
     }
 }

 int Gate::whichChildAdvanced(Node* node, int* ParentIDLR){
     int index=this->whichChild(node);
     if(index==-1){
         for(int i=0; i<this->countChildren(); i++){
             if(this->Children[i]->countChildren()!=0){
                return dynamic_cast<Gate*>(this->Children[i])->whichChildAdvanced(node, ParentIDLR);
             }
         }
     }else{
        *ParentIDLR=this->idLR; 
        return index; //which side of G is node
     }
 }

/**
 * @brief Internal function that produces path probabilities for one point
 * 
 * @param node node where the path ends
 * @param X relevant row of the design matrix X (one point)
 * @param result double to store result in
 * @return double path probability value
 */
 double Gate::getPathProb_internal(Node* node, mat X, double result){ //just one i at a time
     //cout<<"Inside Internal"<<endl;
     //cout<<"This idLR: "<<this->idLR<<endl;
     //cout<<"Node idLR: "<<node->idLR<<endl;
     if(node->idLR!=this->idLR){
         //cout<<"The idLR are not the same"<<endl;
         Gate* parent=node->getParent();
         //cout<<"Consider "<<node->name<<endl;
         int   which=parent->whichChild(node);
         //cout<<node->name<<" is the "<<which<<"th child of "<<parent->name<<endl;
         vec pi_helper=vectorise(parent->pi_calculator(X,parent->gamma));
         //pi_helper.print("pi helper:");
         double pi;
         if(which==0){
              pi=as_scalar(1-sum(pi_helper));
         }else{
             pi=as_scalar(pi_helper[which-1]);
         }
         //cout<<"my pi is "<<pi<<endl;
         result=result*pi;
         return this->getPathProb_internal(parent,X,result); //put return here 
     }else{ 
     return result;
     }
 }

/**
 * @brief Function that produces path probabilities for one point
 * 
 * @param node node where the path ends
 * @param X relevant row of the design matrix X (one point)
 * @return double path probability value
 */
 double Gate::getPathProb(Node* node, mat X){ //X is for one i
     //cout<<"Inside External"<<endl;
     double result=1;
     //cout<<"Setting result to be "<<result<<endl;
    return this->getPathProb_internal(node,X,result);
 }

/**
 * @brief Internal function that produces path probabilities for all points
 * 
 * @param node node where the path ends
 * @param X  design matrix X 
 * 
 * @return vec vector containing path probabilities for all points
 */
vec Gate::getPathProb_mat(Node* node, mat X){ //rows are observations and columns are experts
     vec result(X.n_rows);     
     for(int i=0;i<X.n_rows;i++){
         result[i]=this->getPathProb(node,X.row(i));
     }     
     return result;
 }

/**
 * @brief Prediction function for some new data X
 * 
 * @param X new data X
 * @return vec vector of predicted values y hat
 */
 vec Gate::predict(mat X){
     vector<Node*> terminals=this->getTerminalNodes();
     mat helper(X.n_rows,terminals.size());
     for(int i=0;i<terminals.size();i++){
     vec pathprobs=this->getPathProb_mat(terminals[i],X);
     vec est=X*(dynamic_cast<Expert*>(terminals[i])->beta);
     helper.col(i)=pathprobs%est;
     }
     return sum(helper,1);
 }

/**
 * @brief MCMC step which updates gating parameters gamma
 * this function is virtual in Node and will update beta and sigma instead of gamma if called for an expert
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq expert variance value on a log scale (if normal expert)
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 */
 void Gate::MCMC_internal(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final){
        cout<<"Updating gamma for gate "<<name<<endl;
        mat z=this->getZ(z_final);
        mat Omega;
        mat myX=this->subsetX(X,this->getPointIndices(z_final));
        //cout<<"Before: "<<this->gamma<<endl;
        this->gamma=this->updateGamma(this->gamma,myX,z,Omega);
        for(int i=0;i<this->countChildren();i++){
        this->Children[i]->MCMC_internal(y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_final);
        }
        //cout<<"After: "<<this->gamma<<endl;   
 }

/**
 * @brief MCMC step which updates gating parameters gamma and allocations z once 
 * 
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq expert variance value on a log scale (if normal expert)
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 * @return vector<Node*> updated vector of pointers to experts to which each point has been allocated
 */
vector<Node*> Gate::MCMC_OneRun(vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final){
    this->MCMC_internal(y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_final);
    cout<<"Updating allocations"<<endl;
    z_final=this->updateZ(y,X,z_final);
    return z_final;
}

/**
 * @brief MCMC that updates gating parameters gamma, expert parameters gamma and beta, allocations N times
 * 
 * @param N number of times to run the MCMC chain
 * @param y response vector
 * @param X design matrix
 * @param logsigma_sq expert variance value on a log scale (if normal expert)
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 * @return vector<Node*> updated vector of pointers to experts to which each point has been allocated
 */
vector<Node*> Gate::MCMC(int N, vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta, double a, double b, vector<Node*> z_final){
    vector<Node*> z_new=this->MCMC_OneRun(y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_final);
    ofstream f;
    f.open("results.json");
    f << "[";
    for(int i=0;i<N;i++){
        cout<<"Run number "<<i<<endl;
        z_new=this->MCMC_OneRun(y,X,logsigma_sq,mu_beta,Sigma_beta,a,b,z_new);
        f << this->jsonify() << ",";
    }
    //f << "[" << this->jsonify() << ",";
    //f << this->jsonify() << ",";
    f << "]" << endl;
    f.close();
    return z_new;
}

/**
 * @brief Function that tanslates current state of the tree to a string in json format
 * 
 * @param indent spacing variable (always 0 see wrapper below) 
 * @return string describing current state of the tree
 */
string Gate::jsonify(int indent) {
  map<string, string> m;
  string s;
  m["__name__"] = "\""+this->name+"\"";
  m["__type__"] = "\"gate\"";
  m["_gamma_"] = "\n" + mat2arraystring(this->gamma, indent+6);
  s = "\n" + string(indent+4, ' ') + "[";
  string comma = "";
  for (int i=0; i<this->countChildren(); i++) {
      s = s + comma + "\n" + this->Children[i]->jsonify(indent+6);
      comma=",";  
  }
  m["children"] = s + "\n" + string(indent+4, ' ') + "]";
  return jsondict(m, indent);

} 

/**
 * @brief Function that tanslates current state of the tree to a string in json format
 * 
 * @return string describing current state of the tree
 */
string Gate::jsonify() {
  return this->jsonify(0);
} 

/**
 * @brief Updates the allocations z with beta and sigma integrated out (only if expert models are normal models)
 *  
 * @param y response vector
 * @param X design matrix
 * @param z_final vector of length n of pointers to experts to which each point has been allocated
 * @param mu_beta prior mean for beta
 * @param Sigma_beta prior variance-covariance matrix for beta
 * @param a prior shape parameter for Inverse Gamma
 * @param b prior scale parameter for Inverse Gamma 
 * @return vector<Node*> updated vector of pointers to experts to which each point has been allocated
 */
vector<Node*> Gate::updateZparamsIntegratedOut(vec y, mat X, vector<Node*> z_final,vec mu_beta, mat Sigma_beta, double a, double b){
vector<Node*> terminals=this->getTerminalNodes();
vec alpha(terminals.size());
for(int j=0;j<y.size();j++){
    for(int i=0;i<terminals.size();i++){
        Expert* current=dynamic_cast<Expert*>(terminals[i]);
        vector<Node*> z_helper=z_final;
        z_helper[j]=current;
        vec points=current->getPointIndices(z_helper);
        double marginalY=current->expertmodel->logMarginalPosteriorY(this->subsetY(y,points), this->subsetX(X,points), mu_beta, Sigma_beta, a, b);
        double pathProb=log(this->getPathProb(current,X.row(j)));
        alpha[i]=marginalY+pathProb;
}
    z_final[j]=this->updateZ_onepoint_sample(alpha); //alpha is standardised inside there
}
return z_final;
}

/**
 * @brief An internal function which helps to construct a numeric vector describing a tree
 * This function is called at the node level and returns a vector containing the number of children 
 * of the gate. Then, the same is performed for each of the children of the gate. If the child is an expert,
 * the describeTreeInternal function will be called from the expert level
 * @param description vector to be filled with integers describing a tree
 * @return vector<int> pointer to a vector of integers with the number of children of the gate
 */
vector<int> Gate::describeTreeInternal(vector<int>* description){

        description->push_back(static_cast<int>(Children.size()));

        for(int i=0;i<Children.size();i++){
            Children[i]->describeTreeInternal(description);
        }

        return *description;
};

/**
 * @brief Calculates the brute force log likelihood for the whole tree
 * 
 * @param y response vector
 * @param X design matrix
 * @return double log likelihood value
 */
double Gate::totalLogLikelihood(vec y, mat X){
    //cout<<"Starting calculation"<<endl;
    vector<Node*> terminals=this->getTerminalNodes();
    //for(int i=0;i<terminals.size();i++) cout<<terminals[i]->name<<endl;
    int n=static_cast<int>(X.n_rows);
    //cout<<"n "<<n<<endl;
    vec result(n);
    vec result_perpoint(terminals.size());
    for(int j=0;j<n;j++){
        //cout<<"j="<<j<<endl;
        result_perpoint.zeros();
    for(int i=0;i<terminals.size();i++){
        //cout<<"i="<<i<<endl;
        Expert* current=dynamic_cast<Expert*>(terminals[i]);
        //cout<<"Current name "<<current->name<<endl;
        //cout<<"beta "<<current->beta<<endl;
        //cout<<"sigma"<<current->logsigma_sq<<endl;
        //cout<<"idLR"<<current->idLR<<endl;
        vec index(1);
        index.fill(j);
        //mat myX=X.row(j);
        mat myX=this->subsetX(X,index);
        vec myY=this->subsetY(y,index);
        //myX.print("X:");
        //myY.print("y:");
        result_perpoint[i]=this->getPathProb(current,myX)*as_scalar(current->expertmodel->density(myY,myX*(current->beta),current->logsigma_sq));
        //cout<<"result per point "<<result_perpoint[i]<<endl;
    }
        result[j]=log(sum(result_perpoint));
        //result.print("result:");
    }
    //result.print("result");
    //cout<<"Sum "<<sum(result)<<endl;
 return sum(result);  
}

void  Gate::swapMethod(Gate* gate, int which){
    if(this->Parent->Parent==NULL){
        this->swapRoot(gate,which);
    }else{
        int k;
        int m;
        int i;
        int check=this->whichChild(gate);
    if(check==-1){
        i=this->whichSide(gate->Parent); 
        //cout<<"i denotes which child of "<<this->name<<" leads to "<<gate->Parent->name<<endl;
        //cout<<"i: "<<i<<endl;
        k=this->Parent->whichChild(this);
        //cout<<"k denotes which child of "<<this->Parent->name<<" is "<<this->name<<endl;
        //cout<<"k: "<<k<<endl;
        m=gate->Parent->whichChild(gate);
        //cout<<"m denotes which child of "<<gate->Parent->name<<" is "<<gate->name<<endl;
        //cout<<"m: "<<m<<endl;
        //cout<<"Set the k-th child of "<< this->Parent->name<<" to be "<< gate->name<<endl;
        this->Parent->Children[k]=gate;
        //cout<<"Set the m-th child of "<< gate->Parent->name<<" to be "<< this->name<<endl;
        gate->Parent->Children[m]=this;
        //cout<<"Set the parent of "<< gate->name<<" to be "<<this->Parent->name<<endl;
        Gate* temp=gate->Parent; //save before altering
        gate->Parent=this->Parent;
        //cout<<"Set the parent of "<< this->name<<" to be "<<temp->name<<endl;
        this->Parent=temp;
        //cout<<"Set the i-th child of "<<this->name<<" to be "<<gate->Children[which]->name<<endl;
        Gate* temp2=dynamic_cast<Gate*>(this->Children[i]);
        this->Children[i]=gate->Children[which];
        //cout<<"Set the "<<which<<"-th child of "<<gate->name<<" to be "<<temp2->name<<endl;
        gate->Children[which]=temp2;
        //cout<<"Set the parent of "<<gate->Children[which]->name<<" to be "<<gate->name<<endl;
        gate->Children[which]->Parent=gate;
        //cout<<"Set the parent of "<<this->Children[i]->name<<" to be "<<this->name<<endl;
        this->Children[i]->Parent=this;
    }else{
        //cout<<"Dealing with a child and a parent"<<endl;
        k=this->Parent->whichChild(this);
        //cout<<"k denotes which child of "<<this->Parent->name<<" is "<<this->name<<endl;
        //cout<<"k: "<<k<<endl;
        m=gate->Parent->whichChild(gate);
        //cout<<"m denotes which child of "<<gate->Parent->name<<" is "<<gate->name<<endl;
        //cout<<"m: "<<m<<endl;
        //cout<<"Set the k-th child of "<< this->Parent->name<<" to be "<< gate->name<<endl;
        this->Parent->Children[k]=gate;
        //cout<<"Set the "<<which<<"-th child of "<< gate->name<<" to be "<< this->name<<endl;
        //cout<<"Substitute (and save a copy of) "<<gate->Children[which]->name<<" for "<<this->name<<endl;
        Node* helper=gate->Children[which];
        gate->Children[which]=this;
        //cout<<"Set the parent of "<<gate->Children[which]->name<<" to be "<< gate->name<<endl;
        Gate* helper0=gate->Children[which]->Parent;
        gate->Children[which]->Parent=gate;
        //cout<<"Set the parent of "<<gate->name<<" to be "<< helper0->name<<endl;
        gate->Parent=helper0;
        //cout<<"Set the "<<m<<"-th child of "<< this->name <<" to be the omitted "<<helper->name<<endl;
        this->Children[m]=helper;
        //cout<<"Set the parent of "<<this->Children[m]->name<<" to be "<<this->name<<endl;
        this->Children[m]->Parent=this;
        }
    }
}

void Gate::swapRoot(Gate* gate, int which){
        //cout<<"Swapping root node"<<endl;
        int check=this->whichChild(gate);
        int m;
        int i;
        Gate* RootParent=this->Parent;
        if(check==-1){
            i=this->whichSide(gate->Parent); 
            //cout<<"i denotes which child of "<<this->name<<" leads to "<<gate->Parent->name<<endl;
            //cout<<"i: "<<i<<endl;
            m=gate->Parent->whichChild(gate);
            //cout<<"m denotes which child of "<<gate->Parent->name<<" is "<<gate->name<<endl;
            //cout<<"m: "<<m<<endl;
            //cout<<"Set the m-th child of "<< gate->Parent->name<<" to be "<< this->name<<endl;
            gate->Parent->Children[m]=this;
            //cout<<"Set the parent of "<<this->name<<" to be "<<gate->Parent->name<<endl;
            this->Parent=gate->Parent;
            //cout<<"Set the "<<which<<"-th child of "<<gate->name<<" to be "<<this->Parent->name<<endl;
            Node* omitted=gate->Children[which];
            gate->Children[which]=this->Parent;
            //cout<<"Set the "<<i<<"-th child of "<<this->name<<" to be "<<omitted->name<<endl;
            this->Children[i]=omitted;
            //cout<<"Set the parent of "<<gate->Children[which]->name<<" to be "<<gate->name<<endl;
            gate->Children[which]->Parent=gate;
            //cout<<"Set the parent of "<<this->Children[i]->name<<" to be "<<this->name<<endl;
            this->Children[i]->Parent=this;
            //cout<<"Set the parent of "<<gate->name<<" to be NULL"<<endl;
            gate->Parent=RootParent;
            gate->Parent->Children[0]=gate;
        }else{
            //cout<<"Swapping root node and its child gate"<<endl;
            m=check;
            //cout<<"m denotes which child of "<<gate->Parent->name<<" is "<<gate->name<<endl;
            //cout<<"m: "<<m<<endl;
            //cout<<"Set the "<<which<< "-th child of "<< gate->name<<" to be "<< this->name<<endl;
            Node* helper=gate->Children[which];
            gate->Children[which]=this;
            //cout<<"Set the parent of "<<this->name<<" to be "<<gate->name<<endl;
            this->Parent=gate;
            //cout<<"Set the "<<m<<"-th chid of "<<this->name<<" to be "<<helper->name<<endl;
            this->Children[m]=helper;
            //cout<<"Set the parent of "<<this->Children[m]->name<<" to be "<<this->name<<endl;
            this->Children[m]->Parent=this;
            //cout<<"Set the parent of "<<gate->name<<" to be NULL"<<endl;
            gate->Parent=RootParent;
            gate->Parent->Children[0]=gate;
        }
}

void Gate::swap(Gate* gate, int which, vec y, mat X){
    cout<<"Copy from "<<this->name<<" downwards"<<endl;
    Gate* backup=this->copyStructure();
    cout<<"Save the pointer to parent "<<this->Parent->name<<endl;
    Gate* backupParent=this->Parent;
    int   i=this->Parent->whichChild(this);
    //cout<<this->name<<" is the "<<i<<"-th child of "<<this->Parent->name<<endl;
    cout<<"Swaping "<<this->name<<" with "<<gate->name<<endl;
    //cout<<"Most senior gate in this tree is "<<this->mostSeniorGate()->name<<endl;
    double loglik_old=this->mostSeniorGate()->totalLogLikelihood(y,X);
    this->swapMethod(gate,which);
    double loglik_new=this->mostSeniorGate()->totalLogLikelihood(y,X);
    //Might want to add some code which only calculates likelihood for the relevant segment. 
    cout<<"Loglik old "<<loglik_old<<endl;
    cout<<"Loglik new "<<loglik_new<<endl;
    double acceptance=loglik_new-loglik_old;
    double u=randu();
    bool accept=u<exp(acceptance);
    if(accept==1){
       this->mostSeniorGate()->issueID();
       this->mostSeniorGate()->issueIDLR();
       cout<<"Swap has been accepted"<<endl;
    }else{
       cout<<"Swap has been rejected"<<endl;
       cout<<"Want to reinstate the original structure"<<endl;
       cout<<"Set the "<<i<<"-th child of "<<backupParent->name<<" to be "<<backup->name<<" from the backup hard copy"<<endl;
       backupParent->Children[i]=backup;
       cout<<"Set the parent of "<<backupParent->Children[i]->name<<" to be "<<backupParent->name<<endl;
       backupParent->Children[i]->Parent=backupParent;
    }
}

void  Gate::makeThisRootParent(Node* root){
    this->Parent=NULL;
    this->addChild(root);
    root->Parent=this;
}

void  Gate::estimateAllBetas(vector<Node*> z_assign,vec y, mat X, double logsigma_sq, vec mu_beta, mat Sigma_beta){
    vector<Node*> terminals=this->getTerminalNodes();
    for(int i=0;i<terminals.size();i++){
        Expert* current=dynamic_cast<Expert*>(terminals[i]);
        current->beta=current->expertmodel->findBeta(current->subsetY(y,current->getPointIndices(z_assign)),current->subsetX(X,current->getPointIndices(z_assign)),logsigma_sq,mu_beta,Sigma_beta);
    }
}
   
void  Gate::setAllSigmas(double value){
    vector<Node*> terminals=this->getTerminalNodes();
    for(int i=0;i<terminals.size();i++){
        Expert* current=dynamic_cast<Expert*>(terminals[i]);
        current->logsigma_sq=value;
    }
}

void  Gate::estimateAllGamas(vector<Node*> z_assign, mat X, mat Omega){
    mat myz=this->getZ(z_assign);
    this->gamma=this->findGammaMLE(this->subsetX(X,this->getPointIndices(z_assign)),myz,Omega);
    vector<Node*> desc=this->getDescendants();
    for(int i=0;i<desc.size();i++){
        if(desc[i]->countChildren()!=0){
            Gate* current=dynamic_cast<Gate*>(desc[i]);
            mat z=current->getZ(z_assign); 
            current->gamma=current->findGammaMLE(current->subsetX(X,current->getPointIndices(z_assign)),z,Omega);
        }
    }
}