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
#include "NormalFamily.h"

using namespace std;
using namespace arma;

/**
 * @brief Construct a new Gate:: Gate object
 * 
 */
Gate::Gate(){
    Parent=NULL;
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
    //cout<<"Replace the "<<which<<"-th child of "<<this->name<<endl;
    //cout<<"Replace "<< this->Children[which]->name <<" by "<<newChild->name<<endl;
    //cout<<"Set "<<this->name<<" to be the parent of "<<newChild->name<<endl;
    //cout<<"Remove the parent of "<< this->Children[which]->name<<endl;
    this->Children[which]->Parent=NULL;
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

vector<Node*> Gate::getGates(){
    vector<Node*> desc=this->getDescendants();
    vector<Node*> result;
    result.push_back(this);
    for(int i=0; i<desc.size();i++){
        if(desc[i]->countChildren()!=0) result.push_back(desc[i]);
    }
    return result;
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

int Gate::countGates(){
    vector<Node*> desc=this->getDescendants();
    vector<Node*> terminals=this->getTerminalNodes();
    return 1+desc.size()-terminals.size();
}

int Gate::countTerminals(){
    vector<Node*> terminals=this->getTerminalNodes();
    return terminals.size();
}

int Gate::getMaxGateID(){
vector<Node*> gates=this->getGates();
vec gate_ids(gates.size());
for(int i=0;i<gates.size();i++) gate_ids[i]=gates[i]->id;
return static_cast<int>(gate_ids.max());
}


int Gate::getMaxExpertID(){
vector<Node*> experts=this->getTerminalNodes();
vec expert_ids(experts.size());
for(int i=0;i<experts.size();i++) expert_ids[i]=experts[i]->id;
return static_cast<int>(expert_ids.max());
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

double Gate::loglik_complete(vec y, mat X, vector<Node*> z_assign){
   vector<Node*> terminals=this->getTerminalNodes();
   double result=0;
   for(int j=0;j<terminals.size();j++){
       Expert* myExpert=dynamic_cast<Expert*>(terminals[j]);
       vec myY=this->subsetY(y,myExpert->getPointIndices(z_assign));
       mat myX=this->subsetX(X,myExpert->getPointIndices(z_assign));
       for(int i=0;i<myY.size();i++){
            vec y_helper(1);
            y_helper[0]=myY[i];
            result+=sum(log(this->getPathProb(myExpert,myX.row(i)))+myExpert->expertmodel->logdensity(y_helper,myX.row(i)*myExpert->beta,myExpert->logsigma_sq));
        }
   }
    return result;
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
    //vec diagonals(p*r);
    //diagonals.fill(0.00001);
    //Omega=diagmat(diagonals);
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
    //vec diagonals(p*r);
    //diagonals.fill(0.00001);
    //Omega=diagmat(diagonals);
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
   //diagonals.fill(0.00001);
   //Omega=diagmat(diagonals);
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
        //cout<<"Updating point "<<i<<endl;
        vec index(1);
        index.fill(points[i]);
        int current=static_cast<int>(points[i]);
        vec alpha=this->getSampleProbsForZ(this->subsetY(y,index),X.row(current));
        //cout<<alpha<<endl;
        result[current]=this->updateZ_onepoint_sample(alpha);
        //cout<<"Allocated to "<<result[current]->name<<endl;
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
         alpha[i]=this->getPathProb(current,X)*sum(current->expertmodel->density(y,X*beta,logsigma_sq));
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
        //alpha.print("alphas:");
        //cout<<"rnum "<<rnum<<endl;
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
     if(this->idLR==node->idLR){
        //cout<<"I am done "<<endl;
        return result;
     } 
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
         this->getPathProb_internal(parent,X,result); //put return here 
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
     //cout<<dynamic_cast<Expert*>(terminals[i])->beta<<endl;
     //cout<<pathprobs<<endl;
     vec est=X*(dynamic_cast<Expert*>(terminals[i])->beta); //not true for GLM
     helper.col(i)=pathprobs%est;
     }
     //helper.print("helper:");
     //vec final=sum(helper,1);
     //final.print("summed up:");
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
 void Gate::MCMC_internal(vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega, vector<Node*> z_final){
        //cout<<"Updating gamma for gate "<<name<<endl;
        mat z=this->getZ(z_final);
        mat myX=this->subsetX(X,this->getPointIndices(z_final));
        //cout<<"Before: "<<this->gamma<<endl;
        this->gamma=this->updateGamma(this->gamma,myX,z,Omega);
        //cout<<"After: "<<this->gamma<<endl; 
        for(int i=0;i<this->countChildren();i++){
        this->Children[i]->MCMC_internal(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
        }
          
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
vector<Node*> Gate::MCMC_OneRun(vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega, vector<Node*> z_final){
    //cout<<"I am in"<<endl;
    this->MCMC_internal(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
    //cout<<"Updating allocations"<<endl;
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
vector<Node*> Gate::MCMC(int N, vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega, vector<Node*> z_final){
    //cout<<"Welcome to MCMC"<<endl;
    //cout<<"Perform an initial MCMC run"<<endl;
    vector<Node*> z_new=this->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
    ofstream f;
    f.open("results.json");
    f << "[";
    for(int i=0;i<N;i++){
       //cout<<"Run number "<<i<<endl;
       z_new=this->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_new);
       f << this->jsonify() << ",";
    }
    //f << "[" << this->jsonify() << ",";
    //f << this->jsonify() << ",";
    f << "]" << endl;
    f.close();
    return z_new;
}

vector<Node*> Gate::MCMC(int N, vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega, vector<Node*> z_final, mat* Predictions, mat Xnew, mat* gammas1, mat* gammas2){
    vector<Node*> z_new=this->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
    (*gammas1).col(0)=this->gamma;
    (*gammas2).col(0)=dynamic_cast<Gate*>(this->Children[1])->gamma;
    ofstream f;
    f.open("results.json");
    f << "[";
    for(int i=0;i<N;i++){
        //cout<<"Run number "<<i<<endl;
        z_new=this->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_new);
        (*gammas1).col(i+1)=this->gamma;
        (*gammas2).col(i+1)=dynamic_cast<Gate*>(this->Children[1])->gamma;
        f << this->jsonify() << ",";
        mat helper=this->predict(Xnew);
        helper.reshape(helper.n_cols,helper.n_rows);
        //cout<<helper<<endl;
        //cout<<this->predict(Xnew)<<endl;
        (*Predictions).row(i)=helper;
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
     //alpha is standardised inside there
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
        //cout<<"result for point "<< j<<" is "<< result_perpoint<<endl;
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
    accept=1;
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
        //cout<<"Estimating for "<<current->name<<endl;
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

double Gate::dnorm(double y, double mu, double sigma_sq){
     return exp(-0.5*(log(2*M_PI)+log(sigma_sq))-pow(y-mu,2)/(2*sigma_sq));
}
    
void Gate::splitEmptyExpert(Expert* ExpertToSplit, Expert* ExpertToAdd, Gate* GateToAdd, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega, double* prior_betastar, double* prior_sigmastar, double* prior_gamma, double* q_betastar,double* q_sigmastar, double* q_gamma){
    vec mu_zeros(mu_beta.size()); mu_zeros.fill(0);
    //1) Draw all parameters from the priors
    ExpertToAdd->beta=mvnrnd(mu_beta, Sigma_beta,1);
    ExpertToSplit->beta=mvnrnd(mu_beta, Sigma_beta,1);
    ExpertToAdd->logsigma_sq=log(1/randg( distr_param(a,1/b)));
    ExpertToSplit->logsigma_sq=log(1/randg( distr_param(a,1/b)));
    GateToAdd->gamma=mvnrnd(mu_zeros, Omega,1); 
    //2) Evaluate required results
    *prior_betastar=sum(this->logmvndensity(ExpertToAdd->beta,mu_beta,Sigma_beta))+
                       sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    *prior_sigmastar=ExpertToAdd->expertmodel->IG_log(exp(ExpertToAdd->logsigma_sq),a,b)+
                        ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    *prior_gamma=sum(this->logmvndensity(GateToAdd->gamma,mu_zeros,Omega));
    *q_betastar=*prior_betastar;
    *q_sigmastar=*prior_sigmastar;
    *q_gamma=*prior_gamma;
}

vec Gate::proposeGammaSplit(vec y_sub, mat X_sub,vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon){
    //1) Draw one point at random
    int n_rand=rand() % y_sub.size();  
    mat x_star=X_sub.row(n_rand); x_star.shed_col(0);
    //2) Draw a value for gamma_1
    mat gamma_1=mvnrnd(mu_gamma1, Sigma_gamma1, 1);
    //3) Draw a value for epsilon 
    vec mu_e(1); mu_e.fill(0);
    mat sigma_e(1,1); sigma_e.row(0)=sigma_epsilon;
    double epsilon=as_scalar(mvnrnd(mu_e, sigma_e, 1));
    //4) Infer gamma_0
    mat gamma_0(1,1); gamma_0.row(0)=-as_scalar(x_star*gamma_1)+epsilon;
    //5) Stick both into a vector
    mat gamma=gamma_1;
    gamma.insert_rows(0,gamma_0); 

    return gamma;
}

vector<Node*> Gate::proposeZafterSplit(vec y_sub, mat X_sub, vec points, Gate* GateToAdd, vector<Node*> z_assign, vec* q_z_helper){
    vec alpha(GateToAdd->countChildren());

    for(int i=0;i<y_sub.size();i++){
        vec index(1);
        index.fill(points[i]);
        int current=static_cast<int>(points[i]);
        for(int j=0;j<GateToAdd->countChildren();j++){
            alpha[j]=GateToAdd->getPathProb(GateToAdd->Children[j],X_sub.row(i));
        }
        double sums=sum(alpha);
        alpha=alpha/sums; 
        z_assign[current]=GateToAdd->updateZ_onepoint_sample(alpha);
        (*q_z_helper)[i]=alpha[GateToAdd->whichChild(z_assign[current])];
    }
 
    return z_assign;
}

void Gate::proposeExpertParamsSplit(vec y, mat X, Expert* myExpert,vec mu_beta, mat Sigma_beta, double a, double b, double* q_betastar, double* q_sigmastar){
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    if(X.n_rows>2){
        //1) If not empty
        //2) Calculate estimates for beta and sigma
        vec betahat=myExpert->expertmodel->findBetaMLE(y,X);
        double sigmahat=myExpert->expertmodel->findLogSigmaSqMLE(y,X,betahat);
        mat Sigma_prop= (X.t()*X).i()*exp(sigmahat);      
        //3) Draw a value for beta cenntred around betahat
        myExpert->beta=betahat+mvnrnd(mu_zeros, Sigma_prop,1);
        //4) Draw a value for sigma given beta
        myExpert->logsigma_sq=myExpert->expertmodel->updateSigma(y,X,myExpert->beta,a,b,X.n_rows);
        //5) Record the densities
        *q_betastar+=sum(this->logmvndensity(myExpert->beta,betahat,Sigma_prop));
        *q_sigmastar+=myExpert->expertmodel->IG_log(exp(myExpert->logsigma_sq),a+X.n_rows/2,b+sum(pow(y-X*myExpert->beta,2))/2);
    }else{
        //1) If empty
        //2) Draw both beta and sigma from prior
        myExpert->beta=mvnrnd(mu_beta, Sigma_beta,1);
        myExpert->logsigma_sq=log(1/randg( distr_param(a,1/b)));
        *q_betastar+=sum(this->logmvndensity(myExpert->beta,mu_beta,Sigma_beta));
        *q_sigmastar+=myExpert->expertmodel->IG_log(exp(myExpert->logsigma_sq),a,b);
    }

}

vector<Node*> Gate::deepCopyAllocations(vector<Node*> z_assign, Gate* backup){
    vector<Node*> terminals=backup->getTerminalNodes();
    vector<Node*> result(z_assign.size());
    result=z_assign;
    for(int i=0;i<z_assign.size();i++){
        for(int j=0; j<terminals.size();j++){
            if(z_assign[i]->id==terminals[j]->id) result[i]=terminals[j];
        }
    }
    return result;
}

vector<Node*> Gate::split(vec y, mat X, int* accept, Expert* ExpertToSplit, Expert* ExpertToAdd, Gate* GateToAdd, vector<Node*> z_assign,  vec mu_beta, mat Sigma_beta, vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, double a, double b, mat Omega){// double mu_jump, double sigma_jump, vec mu_beta, mat Sigma_beta, vec* x_record, vec* gamma_record){
    //1) Record the indices of points that reached the expert we are trying to split
    vec  points=ExpertToSplit->getPointIndices(z_assign);
    //2) Check if the chosen expert is empty
    bool empty=points.size()<=2;
    //3) Subset y and X accordingly
    vec y_sub=ExpertToSplit->subsetY(y,points);
    mat X_sub=ExpertToSplit->subsetX(X,points);
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    
    //4) Record relevant information before splitting
    double q_beta=ExpertToSplit->expertmodel->qBeta(y_sub,X_sub,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,mu_beta,Sigma_beta);
    double q_sigma=ExpertToSplit->expertmodel->qSigma(y_sub,X_sub,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,a,b);
    double prior_beta=sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    double prior_sigma=ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    //5) If chosen expert isn't empty, record its log-likelihood
    double loglik_old=100000000;
    if(empty==0) loglik_old=ExpertToSplit->expertmodel->loglik(y_sub,X_sub*(ExpertToSplit->beta),ExpertToSplit->logsigma_sq);
    //if(empty==0) loglik_old=ExpertToSplit->Parent->loglik_complete(y,X,z_assign);
    
    //6) Create necessary backups in case split is rejected
    Gate* backup=ExpertToSplit->Parent->copyStructure();
    Gate* backupParent=ExpertToSplit->Parent->Parent;
    int   k=ExpertToSplit->Parent->whichChild(ExpertToSplit);
    int   m=backupParent->whichChild(backup);
    vector<Node*> z_backup=ExpertToSplit->Parent->deepCopyAllocations(z_assign,backup);
   
    //7) Perform the split
    ExpertToSplit->Parent->replaceChild(k,GateToAdd);
    GateToAdd->addChild(ExpertToSplit);
    GateToAdd->addChild(ExpertToAdd); 
    ExpertToSplit->mostSeniorGate()->issueIDLR();
    ExpertToSplit->mostSeniorGate()->issueID();

    //8) Set up place holders for post split results
    double acceptance;
    double prior_betastar;
    double prior_sigmastar;
    double prior_gamma;
    double q_betastar=0;
    double q_sigmastar=0;
    double q_gamma;
    double loglik_new;
    vector<Node*> z_new(z_assign.size());

    //9)Make parameter proposals
    if(empty==1){
        //cout<<"Splitting an empty expert"<<endl;
        this->splitEmptyExpert(ExpertToSplit,ExpertToAdd, GateToAdd, mu_beta,Sigma_beta,a,b,Omega,&prior_betastar,&prior_sigmastar,&prior_gamma,&q_betastar,&q_sigmastar,&q_gamma);
        acceptance=prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma-
                   prior_beta-prior_sigma-q_betastar-q_sigmastar-q_gamma;
        z_new=z_backup;
    }else{
        //10) Propose a gamma
        GateToAdd->gamma=GateToAdd->proposeGammaSplit(y_sub,X_sub,mu_gamma1,Sigma_gamma1,sigma_epsilon);
        //11) Propose allocations after the split
        vec q_z_helper(y_sub.size());
        z_new=this->proposeZafterSplit(y_sub,X_sub, points, GateToAdd,z_assign,&q_z_helper);
        //12) Subset accordingly
        vec y1=ExpertToAdd->subsetY(y,ExpertToAdd->getPointIndices(z_new)); mat X1=ExpertToAdd->subsetX(X,ExpertToAdd->getPointIndices(z_new)); int n1=static_cast<int>(X1.n_rows);
        vec y2=ExpertToSplit->subsetY(y,ExpertToSplit->getPointIndices(z_new)); mat X2=ExpertToSplit->subsetX(X,ExpertToSplit->getPointIndices(z_new));int n2=static_cast<int>(X2.n_rows);
        //13) Propose expert parameters after the split
        this->proposeExpertParamsSplit(y1,X1,ExpertToAdd,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
        this->proposeExpertParamsSplit(y2,X2,ExpertToSplit,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
        //14) Obtain the probability of arriving to the specific z
        double q_z=sum(log(q_z_helper));
        //15) Obtain the density of chosen gamma
        q_gamma=GateToAdd->q_gammaSplit(X_sub,mu_gamma1,Sigma_gamma1,sigma_epsilon);
        //16) Obtain prior densities of the expert parameters
        prior_betastar=sum(this->logmvndensity(ExpertToAdd->beta,mu_beta,Sigma_beta))+
                       sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
        prior_sigmastar=ExpertToAdd->expertmodel->IG_log(exp(ExpertToAdd->logsigma_sq),a,b)+
                        ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
        prior_gamma=sum(this->logmvndensity(GateToAdd->gamma,mu_zeros,Omega));
        //17) Record the new likelihood
        loglik_new=GateToAdd->loglik_complete(y,X,z_new);
        //18) Calculate acceptance
        acceptance=loglik_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma-
                   loglik_old-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_z-q_gamma; 
        //cout<<"q_z"<<q_z<<endl;
    }

    //UNCOMMENT BELLOW IF WANT TO ZOOM IN ON VALUES
    // cout<<"acceptance"<<acceptance<<endl;
    // cout<<"loglik_old"<<loglik_old<<endl;
    // cout<<"loglik_new"<<loglik_new<<endl;
    // cout<<"prior_betastar"<<prior_betastar<<endl;
    // cout<<"prior_sigmastar"<<prior_sigmastar<<endl;
    // cout<<"prior_gamma"<<prior_gamma<<endl;
    // cout<<"q_betastar"<<q_betastar<<endl;
    // cout<<"q_sigmastar"<<q_sigmastar<<endl;
    // cout<<"q_gamma"<<q_gamma<<endl;

    double u=randu();
    bool acc=log(u)<=acceptance;
    if(acc==1){
        //cout<<"Split has been accepted"<<endl;
        (*accept)=1;
        return z_new;
    }else{
        //cout<<"Split has been rejected"<<endl;
        backupParent->Children[m]=backup;
        backupParent->Children[m]->Parent=backupParent;
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueID();
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueIDLR();
        (*accept)=0;
        return z_backup;
    } 
}


vector<Node*> Gate::split(vec y, mat X, Expert* ExpertToSplit, Expert* ExpertToAdd, Gate* GateToAdd, vector<Node*> z_assign,  vec mu_beta, mat Sigma_beta, vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, double a, double b, mat Omega){// double mu_jump, double sigma_jump, vec mu_beta, mat Sigma_beta, vec* x_record, vec* gamma_record){
    //1) Record the indices of points that reached the expert we are trying to split
    vec  points=ExpertToSplit->getPointIndices(z_assign);
    //2) Check if the chosen expert is empty
    bool empty=points.size()<=2;
    //3) Subset y and X accordingly
    vec y_sub=ExpertToSplit->subsetY(y,points);
    mat X_sub=ExpertToSplit->subsetX(X,points);
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    
    //4) Record relevant information before splitting
    double q_beta=ExpertToSplit->expertmodel->qBeta(y_sub,X_sub,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,mu_beta,Sigma_beta);
    double q_sigma=ExpertToSplit->expertmodel->qSigma(y_sub,X_sub,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,a,b);
    double prior_beta=sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    double prior_sigma=ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    //5) If chosen expert isn't empty, record its log-likelihood
    double loglik_old=100000000;
    if(empty==0) loglik_old=ExpertToSplit->expertmodel->loglik(y_sub,X_sub*(ExpertToSplit->beta),ExpertToSplit->logsigma_sq);
    //if(empty==0) loglik_old=ExpertToSplit->Parent->loglik_complete(y,X,z_assign);
    
    //6) Create necessary backups in case split is rejected
    Gate* backup=ExpertToSplit->Parent->copyStructure();
    Gate* backupParent=ExpertToSplit->Parent->Parent;
    int   k=ExpertToSplit->Parent->whichChild(ExpertToSplit);
    int   m=backupParent->whichChild(backup);
    vector<Node*> z_backup=ExpertToSplit->Parent->deepCopyAllocations(z_assign,backup);
   
    //7) Perform the split
    ExpertToSplit->Parent->replaceChild(k,GateToAdd);
    GateToAdd->addChild(ExpertToSplit);
    GateToAdd->addChild(ExpertToAdd); 
    ExpertToSplit->mostSeniorGate()->issueIDLR();
    ExpertToSplit->mostSeniorGate()->issueID();

    //8) Set up place holders for post split results
    double acceptance;
    double prior_betastar;
    double prior_sigmastar;
    double prior_gamma;
    double q_betastar=0;
    double q_sigmastar=0;
    double q_gamma;
    double loglik_new;
    vector<Node*> z_new(z_assign.size());

    //9)Make parameter proposals
    if(empty==1){
        //cout<<"Splitting an empty expert"<<endl;
        this->splitEmptyExpert(ExpertToSplit,ExpertToAdd, GateToAdd, mu_beta,Sigma_beta,a,b,Omega,&prior_betastar,&prior_sigmastar,&prior_gamma,&q_betastar,&q_sigmastar,&q_gamma);
        acceptance=prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma-
                   prior_beta-prior_sigma-q_betastar-q_sigmastar-q_gamma;
        z_new=z_backup;
    }else{
        //10) Propose a gamma
        GateToAdd->gamma=GateToAdd->proposeGammaSplit(y_sub,X_sub,mu_gamma1,Sigma_gamma1,sigma_epsilon);
        //11) Propose allocations after the split
        vec q_z_helper(y_sub.size());
        z_new=this->proposeZafterSplit(y_sub,X_sub, points, GateToAdd,z_assign,&q_z_helper);
        //12) Subset accordingly
        vec y1=ExpertToAdd->subsetY(y,ExpertToAdd->getPointIndices(z_new)); mat X1=ExpertToAdd->subsetX(X,ExpertToAdd->getPointIndices(z_new)); int n1=static_cast<int>(X1.n_rows);
        vec y2=ExpertToSplit->subsetY(y,ExpertToSplit->getPointIndices(z_new)); mat X2=ExpertToSplit->subsetX(X,ExpertToSplit->getPointIndices(z_new));int n2=static_cast<int>(X2.n_rows);
        //13) Propose expert parameters after the split
        this->proposeExpertParamsSplit(y1,X1,ExpertToAdd,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
        this->proposeExpertParamsSplit(y2,X2,ExpertToSplit,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
        //14) Obtain the probability of arriving to the specific z
        double q_z=sum(log(q_z_helper));
        //15) Obtain the density of chosen gamma
        q_gamma=GateToAdd->q_gammaSplit(X_sub,mu_gamma1,Sigma_gamma1,sigma_epsilon);
        //16) Obtain prior densities of the expert parameters
        prior_betastar=sum(this->logmvndensity(ExpertToAdd->beta,mu_beta,Sigma_beta))+
                       sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
        prior_sigmastar=ExpertToAdd->expertmodel->IG_log(exp(ExpertToAdd->logsigma_sq),a,b)+
                        ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
        prior_gamma=sum(this->logmvndensity(GateToAdd->gamma,mu_zeros,Omega));
        //17) Record the new likelihood
        loglik_new=GateToAdd->loglik_complete(y,X,z_new);
        //18) Calculate acceptance
        acceptance=loglik_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma-
                   loglik_old-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_z-q_gamma; 
        //cout<<"q_z"<<q_z<<endl;
    }

    //UNCOMMENT BELLOW IF WANT TO ZOOM IN ON VALUES
    // cout<<"acceptance"<<acceptance<<endl;
    // cout<<"loglik_old"<<loglik_old<<endl;
    // cout<<"loglik_new"<<loglik_new<<endl;
    // cout<<"prior_betastar"<<prior_betastar<<endl;
    // cout<<"prior_sigmastar"<<prior_sigmastar<<endl;
    // cout<<"prior_gamma"<<prior_gamma<<endl;
    // cout<<"q_betastar"<<q_betastar<<endl;
    // cout<<"q_sigmastar"<<q_sigmastar<<endl;
    // cout<<"q_gamma"<<q_gamma<<endl;

    double u=randu();
    bool acc=log(u)<=acceptance;
    if(acc==1){
        //cout<<"Split has been accepted"<<endl;
        return z_new;
    }else{
        //cout<<"Split has been rejected"<<endl;
        backupParent->Children[m]=backup;
        backupParent->Children[m]->Parent=backupParent;
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueID();
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueIDLR();
        return z_backup;
    } 
}

vector<Node*> Gate::split_at_root(vec y, mat X, int* accept, Expert* ExpertToSplit, Expert* ExpertToAdd, vector<Node*> z_assign,  vec mu_beta, mat Sigma_beta, vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, double a, double b, mat Omega){// double mu_jump, double sigma_jump, vec mu_beta, mat Sigma_beta, vec* x_record, vec* gamma_record){
    //1) Record the indices of all points as integers
    vec  points=ExpertToSplit->getPointIndices(z_assign);
    //2) Set up helpers
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    //3) Record relevant information before splitting
    double q_beta=ExpertToSplit->expertmodel->qBeta(y,X,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,mu_beta,Sigma_beta);
    double q_sigma=ExpertToSplit->expertmodel->qSigma(y,X,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,a,b);
    double prior_beta=sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    double prior_sigma=ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    //4) Record the log-likelihood
    double loglik_old=ExpertToSplit->expertmodel->loglik(y,X*(ExpertToSplit->beta),ExpertToSplit->logsigma_sq);
    //5) Create necessary backups in case split is rejected
    Gate* backup=ExpertToSplit->Parent->copyStructure(); //this is the root gate
    Gate* backupParent=ExpertToSplit->Parent->Parent; //this is the root parent
    vector<Node*> z_backup=ExpertToSplit->Parent->deepCopyAllocations(z_assign,backup);
    //6) Perform the split
    ExpertToSplit->Parent->addChild(ExpertToAdd); //add the second expert as a child of root gate
    ExpertToSplit->mostSeniorGate()->issueIDLR();
    ExpertToSplit->mostSeniorGate()->issueID();
    //7) Set up place holders for post split results
    double acceptance;
    double prior_betastar;
    double prior_sigmastar;
    double prior_gamma;
    double q_betastar=0;
    double q_sigmastar=0;
    double q_gamma;
    double loglik_new;
    vector<Node*> z_new(z_assign.size());

    //8) Propose a gamma
    ExpertToSplit->Parent->gamma=ExpertToSplit->Parent->proposeGammaSplit(y,X,mu_gamma1,Sigma_gamma1,sigma_epsilon);
    //9) Propose allocations after the split
    vec q_z_helper(y.size());
    z_new=this->proposeZafterSplit(y,X, points,ExpertToSplit->Parent,z_assign,&q_z_helper);
    //10) Subset accordingly
    vec y1=ExpertToAdd->subsetY(y,ExpertToAdd->getPointIndices(z_new)); mat X1=ExpertToAdd->subsetX(X,ExpertToAdd->getPointIndices(z_new)); int n1=static_cast<int>(X1.n_rows);
    vec y2=ExpertToSplit->subsetY(y,ExpertToSplit->getPointIndices(z_new)); mat X2=ExpertToSplit->subsetX(X,ExpertToSplit->getPointIndices(z_new));int n2=static_cast<int>(X2.n_rows);
    //11) Propose expert parameters after the split
    this->proposeExpertParamsSplit(y1,X1,ExpertToAdd,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    this->proposeExpertParamsSplit(y2,X2,ExpertToSplit,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    //12) Obtain the probability of arriving to the specific z
    double q_z=sum(log(q_z_helper));
    //13) Obtain the density of chosen gamma
    q_gamma=ExpertToSplit->Parent->q_gammaSplit(X,mu_gamma1,Sigma_gamma1,sigma_epsilon);
    //14) Obtain prior densities of the expert parameters
    prior_betastar=sum(this->logmvndensity(ExpertToAdd->beta,mu_beta,Sigma_beta))+
                   sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    prior_sigmastar=ExpertToAdd->expertmodel->IG_log(exp(ExpertToAdd->logsigma_sq),a,b)+
                    ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    prior_gamma=sum(this->logmvndensity(ExpertToSplit->Parent->gamma,mu_zeros,Omega));
    //15) Record the new likelihood
    loglik_new=ExpertToSplit->Parent->loglik_complete(y,X,z_new);
    //16) Calculate acceptance
    acceptance=loglik_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma-
               loglik_old-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_z-q_gamma; 
    //cout<<"q_z"<<q_z<<endl;
    
    //UNCOMMENT BELLOW IF WANT TO ZOOM IN ON VALUES
    // cout<<"acceptance"<<acceptance<<endl;
    // cout<<"loglik_old"<<loglik_old<<endl;
    // cout<<"loglik_new"<<loglik_new<<endl;
    // cout<<"prior_betastar"<<prior_betastar<<endl;
    // cout<<"prior_sigmastar"<<prior_sigmastar<<endl;
    // cout<<"prior_gamma"<<prior_gamma<<endl;
    // cout<<"q_betastar"<<q_betastar<<endl;
    // cout<<"q_sigmastar"<<q_sigmastar<<endl;
    // cout<<"q_gamma"<<q_gamma<<endl;

    double u=randu();
    bool acc=log(u)<=acceptance;
    if(acc==1){
        //cout<<"Split has been accepted"<<endl;
        (*accept)=1;
        return z_new;
    }else{
        //cout<<"Split has been rejected"<<endl;
        backupParent->Children[0]=backup;
        backupParent->Children[0]->Parent=backupParent;
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueID();
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueIDLR();
        (*accept)=0;
        return z_backup;
    } 
}

vector<Node*> Gate::split_at_root(vec y, mat X, Expert* ExpertToSplit, Expert* ExpertToAdd, vector<Node*> z_assign,  vec mu_beta, mat Sigma_beta, vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, double a, double b, mat Omega){// double mu_jump, double sigma_jump, vec mu_beta, mat Sigma_beta, vec* x_record, vec* gamma_record){
    //1) Record the indices of all points as integers
    vec  points=ExpertToSplit->getPointIndices(z_assign);
    //2) Set up helpers
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    //3) Record relevant information before splitting
    double q_beta=ExpertToSplit->expertmodel->qBeta(y,X,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,mu_beta,Sigma_beta);
    double q_sigma=ExpertToSplit->expertmodel->qSigma(y,X,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,a,b);
    double prior_beta=sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    double prior_sigma=ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    //4) Record the log-likelihood
    double loglik_old=ExpertToSplit->expertmodel->loglik(y,X*(ExpertToSplit->beta),ExpertToSplit->logsigma_sq);
    //5) Create necessary backups in case split is rejected
    Gate* backup=ExpertToSplit->Parent->copyStructure(); //this is the root gate
    Gate* backupParent=ExpertToSplit->Parent->Parent; //this is the root parent
    vector<Node*> z_backup=ExpertToSplit->Parent->deepCopyAllocations(z_assign,backup);
    //6) Perform the split
    ExpertToSplit->Parent->addChild(ExpertToAdd); //add the second expert as a child of root gate
    ExpertToSplit->mostSeniorGate()->issueIDLR();
    ExpertToSplit->mostSeniorGate()->issueID();
    //7) Set up place holders for post split results
    double acceptance;
    double prior_betastar;
    double prior_sigmastar;
    double prior_gamma;
    double q_betastar=0;
    double q_sigmastar=0;
    double q_gamma;
    double loglik_new;
    vector<Node*> z_new(z_assign.size());

    //8) Propose a gamma
    ExpertToSplit->Parent->gamma=ExpertToSplit->Parent->proposeGammaSplit(y,X,mu_gamma1,Sigma_gamma1,sigma_epsilon);
    //9) Propose allocations after the split
    vec q_z_helper(y.size());
    z_new=this->proposeZafterSplit(y,X, points,ExpertToSplit->Parent,z_assign,&q_z_helper);
    //10) Subset accordingly
    vec y1=ExpertToAdd->subsetY(y,ExpertToAdd->getPointIndices(z_new)); mat X1=ExpertToAdd->subsetX(X,ExpertToAdd->getPointIndices(z_new)); int n1=static_cast<int>(X1.n_rows);
    vec y2=ExpertToSplit->subsetY(y,ExpertToSplit->getPointIndices(z_new)); mat X2=ExpertToSplit->subsetX(X,ExpertToSplit->getPointIndices(z_new));int n2=static_cast<int>(X2.n_rows);
    //11) Propose expert parameters after the split
    this->proposeExpertParamsSplit(y1,X1,ExpertToAdd,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    this->proposeExpertParamsSplit(y2,X2,ExpertToSplit,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    //12) Obtain the probability of arriving to the specific z
    double q_z=sum(log(q_z_helper));
    //13) Obtain the density of chosen gamma
    q_gamma=ExpertToSplit->Parent->q_gammaSplit(X,mu_gamma1,Sigma_gamma1,sigma_epsilon);
    //14) Obtain prior densities of the expert parameters
    prior_betastar=sum(this->logmvndensity(ExpertToAdd->beta,mu_beta,Sigma_beta))+
                   sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    prior_sigmastar=ExpertToAdd->expertmodel->IG_log(exp(ExpertToAdd->logsigma_sq),a,b)+
                    ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    prior_gamma=sum(this->logmvndensity(ExpertToSplit->Parent->gamma,mu_zeros,Omega));
    //15) Record the new likelihood
    loglik_new=ExpertToSplit->Parent->loglik_complete(y,X,z_new);
    //16) Calculate acceptance
    acceptance=loglik_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma-
               loglik_old-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_z-q_gamma; 
    //cout<<"q_z"<<q_z<<endl;
    
    //UNCOMMENT BELLOW IF WANT TO ZOOM IN ON VALUES
    // cout<<"acceptance"<<acceptance<<endl;
    // cout<<"loglik_old"<<loglik_old<<endl;
    // cout<<"loglik_new"<<loglik_new<<endl;
    // cout<<"prior_betastar"<<prior_betastar<<endl;
    // cout<<"prior_sigmastar"<<prior_sigmastar<<endl;
    // cout<<"prior_gamma"<<prior_gamma<<endl;
    // cout<<"q_betastar"<<q_betastar<<endl;
    // cout<<"q_sigmastar"<<q_sigmastar<<endl;
    // cout<<"q_gamma"<<q_gamma<<endl;

    double u=randu();
    bool acc=log(u)<=acceptance;
    if(acc==1){
        //cout<<"Split has been accepted"<<endl;
        return z_new;
    }else{
        //cout<<"Split has been rejected"<<endl;
        backupParent->Children[0]=backup;
        backupParent->Children[0]->Parent=backupParent;
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueID();
        dynamic_cast<Gate*>(backupParent->mostSeniorGate())->issueIDLR();
        return z_backup;
    } 
}


double Gate::q_gammaSplit(mat X_sub, vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon){
    //This is the density of proposed gamma in the split. See Architecture Selection chapter for details
    double q_gamma=0;
    vec gamma=this->gamma;
        for(int i=0;i<X_sub.n_rows;i++){
            mat mean_helper(gamma.size(),1);
            mean_helper=mu_gamma1;
            mat x_helper=X_sub.row(i);
            x_helper.shed_col(0);
            mean_helper.insert_rows(0,-x_helper*mu_gamma1);
            mat Sigma11=sigma_epsilon+x_helper*Sigma_gamma1*x_helper.t();
            mat Sigma12=-x_helper*Sigma_gamma1;
            mat Sigma_helper=join_cols(join_rows(Sigma11,Sigma12),join_rows(Sigma12.t(),Sigma_gamma1));
            q_gamma+=sum(exp(this->logmvndensity(gamma,mean_helper,Sigma_helper)));
        }
        q_gamma=log(q_gamma/X_sub.n_rows);
        return q_gamma;
}


vector<Node*> Gate::merge(vec y, mat X, int* accept, Gate* GateToMerge, vector<Node*> z_assign, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega){
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    //1) Determine which experts to merge
    vector<Node*> ExpertsToMerge=GateToMerge->getTerminalNodes();
    //2) Subset the points that have reached the gate
    vec points=GateToMerge->getPointIndices(z_assign);
    //3) Record the log likelihood, priors and parameter densities for the current architecture
    double loglik_old=GateToMerge->loglik_complete(y,X,z_assign); 
    double prior_beta=0;
    double prior_sigma=0;
    double q_beta=0;
    double q_sigma=0;
    double q_gamma=GateToMerge->q_gammaMerge(X,z_assign,Omega); //subsets X based on z_assign inside
    double prior_gamma=sum(this->logmvndensity(GateToMerge->gamma,mu_zeros,Omega));
    for(int i=0;i<ExpertsToMerge.size();i++){
        Expert* current=dynamic_cast<Expert*>(ExpertsToMerge[i]);
        vec points_helper=current->getPointIndices(z_assign);
        vec y_sub=current->subsetY(y,points_helper);
        mat X_sub=current->subsetX(X,points_helper);
        prior_beta=prior_beta+sum(this->logmvndensity(current->beta,mu_beta,Sigma_beta));
        prior_sigma=prior_sigma+current->expertmodel->IG_log(exp(current->logsigma_sq),a,b);
        q_beta=q_beta+current->expertmodel->qBeta(y_sub,X_sub,current->beta,current->logsigma_sq,mu_beta,Sigma_beta);
        q_sigma=q_sigma+current->expertmodel->qSigma(y_sub,X_sub,current->beta,current->logsigma_sq,a,b);
    }
    //4) Create backups in case merge is rejected
    Gate* backup=GateToMerge->copyStructure();
    Gate* backupParent=GateToMerge->Parent;
    int   k=GateToMerge->Parent->whichChild(GateToMerge);
    int   m=backupParent->whichChild(backup);
    vector<Node*> z_backup=GateToMerge->deepCopyAllocations(z_assign,backup);
    //5)Create variables to be filled in
    bool acc;
    vector<Node*> z_new=z_assign;
    //6) Perform the merge
    if(GateToMerge->Parent->Parent==NULL){
        //a) If merging a root gate, only have one expert in the model
        GateToMerge->deleteChildren();
        GateToMerge->addChild(ExpertsToMerge[0]);
    }else{
        GateToMerge->Parent->replaceChild(k,ExpertsToMerge[0]);
    }
    //7) Create a new pointer to the newly formed expert to simplify notation
    Expert* MergedExpert=dynamic_cast<Expert*>(ExpertsToMerge[0]);
    //8) Issue IDs to the newly created architecture
    MergedExpert->mostSeniorGate()->issueID();
    MergedExpert->mostSeniorGate()->issueIDLR();
    //9) Assign all points to one merged expert
    for(int i=0;i<points.size();i++){
         z_new[static_cast<int>(points[i])]=MergedExpert;
    }
    //10) Subset all points that are in the new merged expert
    mat myX=MergedExpert->subsetX(X,MergedExpert->getPointIndices(z_new));
    vec myY=MergedExpert->subsetY(y,MergedExpert->getPointIndices(z_new));
    //11) Propose post merge parameters
    double q_betastar=0; double q_sigmastar=0;
    MergedExpert->Parent->proposeExpertParamsSplit(myY,myX,MergedExpert,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    //12) Record log-likelihood after the merge
    double loglik_new=MergedExpert->expertmodel->loglik(myY,myX*(MergedExpert->beta),MergedExpert->logsigma_sq);
    //13) Evaluate priors
    double prior_betastar=sum(this->logmvndensity(MergedExpert->beta,mu_beta,Sigma_beta));
    double prior_sigmastar=MergedExpert->expertmodel->IG_log(exp(MergedExpert->logsigma_sq),a,b);
    //14) Calculate the acceptance probability
    double acceptance=loglik_new+prior_betastar+prior_sigmastar+q_beta+q_sigma+q_gamma-
                      loglik_old-prior_beta-prior_sigma-prior_gamma-q_betastar-q_sigmastar;
   //15) Check if accepted
    double u=randu();
    acc=log(u)<acceptance;
    if(acc==1){
         //cout<<"Merge has been accepted"<<endl;
         (*accept)=1;
         return z_new;
    }else{
         //cout<<"Merge has been rejected"<<endl;
         backupParent->Children[m]=backup;
         backupParent->Children[m]->Parent=backupParent;
         backupParent->mostSeniorGate()->issueID();
         backupParent->mostSeniorGate()->issueIDLR();
         (*accept)=0;
         return z_backup;
}
}

vector<Node*> Gate::merge(vec y, mat X, Gate* GateToMerge, vector<Node*> z_assign, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega){
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    //1) Determine which experts to merge
    vector<Node*> ExpertsToMerge=GateToMerge->getTerminalNodes();
    //2) Subset the points that have reached the gate
    vec points=GateToMerge->getPointIndices(z_assign);
    //3) Record the log likelihood, priors and parameter densities for the current architecture
    double loglik_old=GateToMerge->loglik_complete(y,X,z_assign); 
    double prior_beta=0;
    double prior_sigma=0;
    double q_beta=0;
    double q_sigma=0;
    double q_gamma=GateToMerge->q_gammaMerge(X,z_assign,Omega); //subsets X based on z_assign inside
    double prior_gamma=sum(this->logmvndensity(GateToMerge->gamma,mu_zeros,Omega));
    for(int i=0;i<ExpertsToMerge.size();i++){
        Expert* current=dynamic_cast<Expert*>(ExpertsToMerge[i]);
        vec points_helper=current->getPointIndices(z_assign);
        vec y_sub=current->subsetY(y,points_helper);
        mat X_sub=current->subsetX(X,points_helper);
        prior_beta=prior_beta+sum(this->logmvndensity(current->beta,mu_beta,Sigma_beta));
        prior_sigma=prior_sigma+current->expertmodel->IG_log(exp(current->logsigma_sq),a,b);
        q_beta=q_beta+current->expertmodel->qBeta(y_sub,X_sub,current->beta,current->logsigma_sq,mu_beta,Sigma_beta);
        q_sigma=q_sigma+current->expertmodel->qSigma(y_sub,X_sub,current->beta,current->logsigma_sq,a,b);
    }
    //4) Create backups in case merge is rejected
    Gate* backup=GateToMerge->copyStructure();
    Gate* backupParent=GateToMerge->Parent;
    int   k=GateToMerge->Parent->whichChild(GateToMerge);
    int   m=backupParent->whichChild(backup);
    vector<Node*> z_backup=GateToMerge->deepCopyAllocations(z_assign,backup);
    //5)Create variables to be filled in
    bool acc;
    vector<Node*> z_new=z_assign;
    //6) Perform the merge
    if(GateToMerge->Parent->Parent==NULL){
        //a) If merging a root gate, only have one expert in the model
        GateToMerge->deleteChildren();
        GateToMerge->addChild(ExpertsToMerge[0]);
    }else{
        GateToMerge->Parent->replaceChild(k,ExpertsToMerge[0]);
    }
    //7) Create a new pointer to the newly formed expert to simplify notation
    Expert* MergedExpert=dynamic_cast<Expert*>(ExpertsToMerge[0]);
    //8) Issue IDs to the newly created architecture
    MergedExpert->mostSeniorGate()->issueID();
    MergedExpert->mostSeniorGate()->issueIDLR();
    //9) Assign all points to one merged expert
    for(int i=0;i<points.size();i++){
         z_new[static_cast<int>(points[i])]=MergedExpert;
    }
    //10) Subset all points that are in the new merged expert
    mat myX=MergedExpert->subsetX(X,MergedExpert->getPointIndices(z_new));
    vec myY=MergedExpert->subsetY(y,MergedExpert->getPointIndices(z_new));
    //11) Propose post merge parameters
    double q_betastar=0; double q_sigmastar=0;
    MergedExpert->Parent->proposeExpertParamsSplit(myY,myX,MergedExpert,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    //12) Record log-likelihood after the merge
    double loglik_new=MergedExpert->expertmodel->loglik(myY,myX*(MergedExpert->beta),MergedExpert->logsigma_sq);
    //13) Evaluate priors
    double prior_betastar=sum(this->logmvndensity(MergedExpert->beta,mu_beta,Sigma_beta));
    double prior_sigmastar=MergedExpert->expertmodel->IG_log(exp(MergedExpert->logsigma_sq),a,b);
    //14) Calculate the acceptance probability
    double acceptance=loglik_new+prior_betastar+prior_sigmastar+q_beta+q_sigma+q_gamma-
                      loglik_old-prior_beta-prior_sigma-prior_gamma-q_betastar-q_sigmastar;
   //15) Check if accepted
    double u=randu();
    acc=log(u)<acceptance;
    if(acc==1){
         //cout<<"Merge has been accepted"<<endl;
         return z_new;
    }else{
         //cout<<"Merge has been rejected"<<endl;
         backupParent->Children[m]=backup;
         backupParent->Children[m]->Parent=backupParent;
         backupParent->mostSeniorGate()->issueID();
         backupParent->mostSeniorGate()->issueIDLR();
         return z_backup;
}
}

 double Gate::q_gammaMerge(mat X, vector<Node*> z_assign, mat Omega){
     mat z=this->getZ(z_assign);
     mat myX=this->subsetX(X,this->getPointIndices(z_assign));
     mat R;
     vec gammahat = this->findGammaQR(myX, z, Omega,&R);
     mat Sigma=(R.t()*R).i();
     return sum(this->logmvndensity(this->gamma,gammahat,Sigma));
 }

vector<Node*> Gate::MCMC_RJ(int N,  bool doRJ, int RJ_every, int L, mat* accept_RJ, vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega,vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, vector<Node*> z_final, mat X_new, mat* predictions, int predict_every){
    //This MCMC records both acceptance and predictions
    if(doRJ==0) cout<<"Welcome to the MCMC"<<endl;
    if(doRJ==1) cout<<"Welcome to the RJ MCMC!"<<endl;
    cout<<"This MCMC wil run "<<N<<" times."<<endl;
    if(doRJ==1) cout<<"There will be "<<L<<" jumps proposed every "<< RJ_every<<" iteration."<<endl;
    Gate* RootParent=this->mostSeniorGate()->Parent;
    Gate* G0=dynamic_cast<Gate*>(RootParent->Children[0]); 
    int jump_no=0;
    int predict_no=0;
    for(int i=0;i<N;i++){
            //Run MCMC, if only one expert, then only update expert parameters
            if(G0->countTerminals()>1){ 
                z_final=G0->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }else{
                dynamic_cast<Expert*>(G0->Children[0])->MCMC_internal(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }
            if(doRJ==1){ //If the RJ is on
                    if(i%RJ_every==0){ //Check if it is the correct iteration
                        for(int j=0;j<L;j++){   
                            G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                            //cout<<"RJ taking place at iteration "<< i<<" for the "<<j+1<<"-th time"<<endl;
                            int RJ_direction= rand() % 2; //Randomly choose the direction of the jump
                            if(G0->mostSeniorGate()->countTerminals()==1) RJ_direction=0; //If there is only one expert propose to split it
                            if(RJ_direction==0){ //in case of a split
                                //cout<<"Performing a split."<<endl;
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                                vector<Node*> terminals=G0->getTerminalNodes();
                                //arma_rng::set_seed_random();
                                int n_rand=rand() % terminals.size(); //chose expert to split at random
                                Expert* ExpertToSplit = dynamic_cast<Expert*>(terminals[n_rand]);
                                //cout<<"Proposing to split "<<ExpertToSplit->name<<endl;
                                Expert* ExpertToAdd= new Expert(); //create a new sibling expert
                                ExpertToAdd->name= "E"+to_string(G0->getMaxExpertID()+1);
                                NormalFamily* NF=new NormalFamily();
                                ExpertToAdd->expertmodel=NF;
                                int accept;
                                if(ExpertToSplit->Parent->countChildren()==1){
                                    z_final=G0->split_at_root(y,X,&accept,ExpertToSplit,ExpertToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
                                }else{
                                    Gate*   GateToAdd=new Gate(); //create a new gate parent
                                    GateToAdd->name="G"+to_string(G0->getMaxGateID()+1);
                                    //G0->printChildren();
                                    z_final=G0->split(y,X,&accept,ExpertToSplit,ExpertToAdd,GateToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
                                }//closes the usual split when more than 1 child
                                //RootParent->Children[0]->printChildren();
                                (*accept_RJ).submat(span(jump_no,jump_no),span(0,0))=1;
                                (*accept_RJ).submat(span(jump_no,jump_no),span(1,1))=accept;
                                jump_no+=1;
                            }else{
                                //cout<<"Performing a merge"<<endl;
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                                vector<Node*> gates=G0->getGates();
                                gates.erase(gates.begin());
                                int n_gates=G0->countGates()-1;
                                Gate* GateToMerge;
                                if(n_gates!=0){ //if there is more than a root gate
                                    //arma_rng::set_seed_random();
                                    int n_rand=rand() % n_gates;
                                    GateToMerge=dynamic_cast<Gate*>(gates[n_rand]); //choose a gate to merge at random
                                }else{
                                    GateToMerge=G0; //merge the root
                                }
                                //cout<<"Propose to merge the children of "<<GateToMerge->name<<endl;
                                //GateToMerge->printChildren();
                                int accept2;
                                z_final=G0->merge(y,X,&accept2,GateToMerge,z_final,mu_beta,Sigma_beta,a,b,Omega);
                                //RootParent->Children[0]->printChildren();
                                (*accept_RJ).submat(span(jump_no,jump_no),span(0,0))=2;
                                (*accept_RJ).submat(span(jump_no,jump_no),span(1,1))=accept2;
                                jump_no+=1;
                            } //closes the merge
                        } //closes for L
                    } //closes if divisible    
                } //closes if doRJ
                    G0=dynamic_cast<Gate*>(RootParent->Children[0]); //record current root node
                    if(i%predict_every==0){
                        (*predictions).col(predict_no)=G0->predict(X_new);
                        predict_no+=1;
                    } //closes predict_every
                } //closes the first for loop
            return z_final;
        } //closes the function


 vector<Node*> Gate::MCMC_RJ(int N,  bool doRJ, int RJ_every, int L, mat* accept_RJ, vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega,vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, vector<Node*> z_final){
    //This MCMC records acceptance, but not predictions
    if(doRJ==0) cout<<"Welcome to the MCMC"<<endl;
    if(doRJ==1) cout<<"Welcome to the RJ MCMC!"<<endl;
    cout<<"This MCMC wil run "<<N<<" times."<<endl;
    if(doRJ==1) cout<<"There will be "<<L<<" jumps proposed every "<< RJ_every<<" iteration."<<endl;
    Gate* RootParent=this->mostSeniorGate()->Parent;
    Gate* G0=dynamic_cast<Gate*>(RootParent->Children[0]); 
    int jump_no=0;
    int predict_no=0;
    for(int i=0;i<N;i++){
            //Run MCMC, if only one expert, then only update expert parameters
            if(G0->countTerminals()>1){ 
                z_final=G0->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }else{
                dynamic_cast<Expert*>(G0->Children[0])->MCMC_internal(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }
            if(doRJ==1){ //If the RJ is on
                    if(i%RJ_every==0){ //Check if it is the correct iteration
                        for(int j=0;j<L;j++){   
                            G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                            //cout<<"RJ taking place at iteration "<< i<<" for the "<<j+1<<"-th time"<<endl;
                            int RJ_direction= rand() % 2; //Randomly choose the direction of the jump
                            if(G0->mostSeniorGate()->countTerminals()==1) RJ_direction=0; //If there is only one expert propose to split it
                            if(RJ_direction==0){ //in case of a split
                                //cout<<"Performing a split."<<endl;
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                                vector<Node*> terminals=G0->getTerminalNodes();
                                //arma_rng::set_seed_random();
                                int n_rand=rand() % terminals.size(); //chose expert to split at random
                                Expert* ExpertToSplit = dynamic_cast<Expert*>(terminals[n_rand]);
                                //cout<<"Proposing to split "<<ExpertToSplit->name<<endl;
                                Expert* ExpertToAdd= new Expert(); //create a new sibling expert
                                ExpertToAdd->name= "E"+to_string(G0->getMaxExpertID()+1);
                                NormalFamily* NF=new NormalFamily();
                                ExpertToAdd->expertmodel=NF;
                                int accept;
                                if(ExpertToSplit->Parent->countChildren()==1){
                                    z_final=G0->split_at_root(y,X,&accept,ExpertToSplit,ExpertToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
                                }else{
                                    Gate*   GateToAdd=new Gate(); //create a new gate parent
                                    GateToAdd->name="G"+to_string(G0->getMaxGateID()+1);
                                    //G0->printChildren();
                                    z_final=G0->split(y,X,&accept,ExpertToSplit,ExpertToAdd,GateToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
                                }//closes the usual split when more than 1 child
                                //RootParent->Children[0]->printChildren();
                                (*accept_RJ).submat(span(jump_no,jump_no),span(0,0))=1;
                                (*accept_RJ).submat(span(jump_no,jump_no),span(1,1))=accept;
                                jump_no+=1;
                            }else{
                                //cout<<"Performing a merge"<<endl;
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                                vector<Node*> gates=G0->getGates();
                                gates.erase(gates.begin());
                                int n_gates=G0->countGates()-1;
                                Gate* GateToMerge;
                                if(n_gates!=0){ //if there is more than a root gate
                                    //arma_rng::set_seed_random();
                                    int n_rand=rand() % n_gates;
                                    GateToMerge=dynamic_cast<Gate*>(gates[n_rand]); //choose a gate to merge at random
                                }else{
                                    GateToMerge=G0; //merge the root
                                }
                                //cout<<"Propose to merge the children of "<<GateToMerge->name<<endl;
                                //GateToMerge->printChildren();
                                int accept2;
                                z_final=G0->merge(y,X,&accept2,GateToMerge,z_final,mu_beta,Sigma_beta,a,b,Omega);
                                //RootParent->Children[0]->printChildren();
                                (*accept_RJ).submat(span(jump_no,jump_no),span(0,0))=2;
                                (*accept_RJ).submat(span(jump_no,jump_no),span(1,1))=accept2;
                                jump_no+=1;
                            } //closes the merge
                        } //closes for L
                    } //closes if divisible    
                } //closes if doRJ
                    G0=dynamic_cast<Gate*>(RootParent->Children[0]); //record current root node
            } //closes the first for loop
            return z_final;
        } //closes the function


vector<Node*> Gate::MCMC_RJ(int N,  bool doRJ, int RJ_every, int L, vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega,vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, vector<Node*> z_final){
    //This MCMC does not record acceptace or predictions
    if(doRJ==0) cout<<"Welcome to the MCMC"<<endl;
    if(doRJ==1) cout<<"Welcome to the RJ MCMC!"<<endl;
    cout<<"This MCMC wil run "<<N<<" times."<<endl;
    if(doRJ==1) cout<<"There will be "<<L<<" jumps proposed every "<< RJ_every<<" iteration."<<endl;
    Gate* RootParent=this->mostSeniorGate()->Parent;
    Gate* G0=dynamic_cast<Gate*>(RootParent->Children[0]); 
    int jump_no=0;
    for(int i=0;i<N;i++){
            //Run MCMC, if only one expert, then only update expert parameters
            if(G0->countTerminals()>1){ 
                z_final=G0->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }else{
                dynamic_cast<Expert*>(G0->Children[0])->MCMC_internal(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }
            if(doRJ==1){ //If the RJ is on
                    if(i%RJ_every==0){ //Check if it is the correct iteration
                        for(int j=0;j<L;j++){   
                            G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                            //cout<<"RJ taking place at iteration "<< i<<" for the "<<j+1<<"-th time"<<endl;
                            int RJ_direction= rand() % 2; //Randomly choose the direction of the jump
                            if(G0->mostSeniorGate()->countTerminals()==1) RJ_direction=0; //If there is only one expert propose to split it
                            if(RJ_direction==0){ //in case of a split
                                //cout<<"Performing a split."<<endl;
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                                vector<Node*> terminals=G0->getTerminalNodes();
                                //arma_rng::set_seed_random();
                                int n_rand=rand() % terminals.size(); //chose expert to split at random
                                Expert* ExpertToSplit = dynamic_cast<Expert*>(terminals[n_rand]);
                                //cout<<"Proposing to split "<<ExpertToSplit->name<<endl;
                                Expert* ExpertToAdd= new Expert(); //create a new sibling expert
                                ExpertToAdd->name= "E"+to_string(G0->getMaxExpertID()+1);
                                NormalFamily* NF=new NormalFamily();
                                ExpertToAdd->expertmodel=NF;
                                if(ExpertToSplit->Parent->countChildren()==1){
                                    z_final=G0->split_at_root(y,X,ExpertToSplit,ExpertToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
                                }else{
                                    Gate*   GateToAdd=new Gate(); //create a new gate parent
                                    GateToAdd->name="G"+to_string(G0->getMaxGateID()+1);
                                    //G0->printChildren();
                                    z_final=G0->split(y,X,ExpertToSplit,ExpertToAdd,GateToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega);
                                }//closes the usual split when more than 1 child
                                //RootParent->Children[0]->printChildren();
                                jump_no+=1;
                            }else{
                                //cout<<"Performing a merge"<<endl;
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                                vector<Node*> gates=G0->getGates();
                                gates.erase(gates.begin());
                                int n_gates=G0->countGates()-1;
                                Gate* GateToMerge;
                                if(n_gates!=0){ //if there is more than a root gate
                                    //arma_rng::set_seed_random();
                                    int n_rand=rand() % n_gates;
                                    GateToMerge=dynamic_cast<Gate*>(gates[n_rand]); //choose a gate to merge at random
                                }else{
                                    GateToMerge=G0; //merge the root
                                }
                                //cout<<"Propose to merge the children of "<<GateToMerge->name<<endl;
                                //GateToMerge->printChildren();
                                z_final=G0->merge(y,X,GateToMerge,z_final,mu_beta,Sigma_beta,a,b,Omega);
                                //RootParent->Children[0]->printChildren();
                                jump_no+=1;
                            } //closes the merge
                        } //closes for L
                    } //closes if divisible    
                } //closes if doRJ
                    G0=dynamic_cast<Gate*>(RootParent->Children[0]); //record current root node
                } //closes the first for loop
            return z_final;
        } //closes the function
    

 vector<Expert*> Gate::whichEmpty(vector<Node*> z_assign){
    vector<Expert*> result;
    vector<Node*> terminals=this->getTerminalNodes();
    for(int i=0; i<terminals.size();i++){
        vec points=terminals[i]->getPointIndices(z_assign);
        if(points.size()<=2) result.push_back(dynamic_cast<Expert*>(terminals[i]));
    }  
    return result;
 }

 int Gate::areAnyExpEmpty(vector<Node*> z_assign){
    vector<Node*> terminals=this->getTerminalNodes();
    vec result(terminals.size());
    result.fill(0);
    for(int i=0;i<terminals.size();i++){
        vec points=terminals[i]->getPointIndices(z_assign);
        if(points.size()<=2) result[i]=1;
    }
    if(sum(result)==0){
        return 0;
    }else{
        return 1;
    }

 }

mat Gate::extractAllParams(){
    int n_exp=this->mostSeniorGate()->countTerminals();
    int n_gates=this->mostSeniorGate()->countGates();
    vector<Node*> terminals=this->mostSeniorGate()->getTerminalNodes();
    vector<Node*> gates=this->mostSeniorGate()->getGates();
    int p=(dynamic_cast<Gate*>(gates[0])->gamma).size();
    mat result(p+1,n_exp+n_gates);
    result.fill(0);
    for(int i=0; i<terminals.size();i++){
        Expert* current=dynamic_cast<Expert*>(terminals[i]);
        result.submat(span(0,p-1),span(i,i))=current->beta;
        result.submat(span(p,p),span(i,i))=current->logsigma_sq;
    }
    for(int j=0;j<gates.size();j++){
        Gate* current=dynamic_cast<Gate*>(gates[j]);
        result.submat(span(0,p-1),span(n_exp+j,n_exp+j))=current->gamma;
    }
    //result.print("result");
    return result;
}