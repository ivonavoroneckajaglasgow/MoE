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

vector<Node*> Gate::getLastGates(){
    vector<Node*> gates=this->getGates();
    vector<Node*> result;
    for(int i=0;i<gates.size();i++){
        //cout<<"Considering gate "<<gates[i]->name<<endl;
        vector<Node*> temp=gates[i]->getChildren();
        //gates[i]->printChildren();
        int children=0;
        for(int j=0;j<temp.size();j++){
            //cout<<temp[j]->name<<" has "<<temp[j]->countChildren()<<" children "<<endl;
            children+=temp[j]->countChildren();
        }
        if(children==0){
            //cout<<"All children are experts"<<endl;
            result.push_back(gates[i]);
        } 
    }
    //for(int k=0;k<result.size();k++) cout<<result[k]->name<<endl;  
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

int Gate::countEmptyTerminals(vector<Node*> z_assign){
    vector<Node*> terminals=this->getTerminalNodes();
    int result=0;
    for(int i=0;i<terminals.size();i++){
        if(terminals[i]->countPoints(z_assign)<=2) result+=1;
    }
    return result;
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
       //cout<<"beta: "<<myExpert->beta<<endl;
       //cout<<"logsigma: "<<myExpert->logsigma_sq<<endl;
       vec myY=this->subsetY(y,myExpert->getPointIndices(z_assign));
       mat myX=this->subsetX(X,myExpert->getPointIndices(z_assign));
       for(int i=0;i<myY.size();i++){
            vec y_helper(1);
            y_helper[0]=myY[i];
            //cout<<"E"<<myExpert->id<<" for point "<<y[i]<<" path prob is: "<<this->getPathProb(myExpert,myX.row(i))<<endl;
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
    mat eta=X*gamma2;
    mat eta_pos=eta.elem(find(eta > 0));
    mat helper=exp(eta);
    helper.elem(find(eta > 0)).ones();
    mat rowsums=this->getRowSumsMat(helper);
    rowsums.elem(find(eta > 0))=this->getRowSumsMat(exp(-eta_pos));
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
     for(int i=0; i<1000; i++){
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
    for(int i=0; i<1000; i++){
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
    Omega=Omega.i();
    //vec diagonals(p*r);
    //diagonals.fill(0.00001);
    //Omega=diagmat(diagonals);
    for(int i=0; i<1000; i++){
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
    mat Omega_chol=chol(Omega.i());
    vec gamma(p*r);
    gamma.zeros();
    if(n==0){
        *R=Omega_chol;
        return gamma;
    }
    for(int i=0; i<1000; i++){
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
   vec mu_gamma(gammaold.size()); mu_gamma.zeros();
   mat R;
   vec gammahat = this->findGammaQR(X, z, Omega,&R);//The R correspods to variance, there is no inverse of Omega taken
   vec v(gammaold.size(),fill::randn);
   mat RHS(R.n_rows,1); RHS.zeros(); RHS.rows(0,gammahat.size()-1)=v;
   //vec gammanew=gammahat+sqrt(SigmaMultiple)*solve(R,RHS);//random mvn
   //vec gammanew=gammahat+solve(R,RHS);//random mvn
   vec gammanew=gammahat+solve(R,RHS);//random mvn
   double loglik_old=this->loglik(z,this->pi_calculator(X,gammaold));
   double loglik_new=this->loglik(z,this->pi_calculator(X,gammanew));
   //OLD:
   //double proposal_old=sum(this->logmvndensity(gammaold,gammahat,&R));
   //double proposal_new=sum(this->logmvndensity(gammanew,gammahat,&R));
   //try this:
   mat Sigma=(R.t()*R).i();
   double proposal_old=sum(this->logmvndensity(gammaold,gammahat,Sigma));
   double proposal_new=sum(this->logmvndensity(gammanew,gammahat,Sigma));
   double prior_old=sum(this->logmvndensity(gammaold,mu_gamma, Omega));
   double prior_new=sum(this->logmvndensity(gammanew,mu_gamma, Omega));
   double acceptance=loglik_new-loglik_old+proposal_old-proposal_new+prior_new-prior_old;
   double u=randu();
   bool accept=u<exp(acceptance);
    //gammahat.print("gammahat:");
    //gammaold.print("gammaold:");
    //gammanew.print("gammanew:");
  if(accept==1){
     //cout<<"accepted"<<endl;
    return gammanew;
  }else{
     //cout<<"rejected"<<endl;
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
   //cout<<"I am in the one without the pointer"<<endl;
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
//int k=static_cast<int>((*R).n_rows);
//return -k/2*log(2*M_PI)+0.5*sum(log(pow((*R).diag(),2)))-0.5*(response-mean).t()*((*R).t()*(*R))*(response-mean);
//cout<<"I am in the one with the pointer"<<endl;
mat result;
mat Sigma=((*R).t()*(*R)).i();
//result =-k/2*log(2*M_PI)+sum(log((*R).diag()))-0.5*sum(pow((*R)*(response-mean),2));
result=this->logmvndensity(response,mean,Sigma);
return vectorise(result);
}

double Gate::poissondensity(int no_expt, double lambda){
    return pow(lambda,no_expt)/ tgamma(no_expt+1)*exp(-lambda);
}

double Gate::geomdensity(int no_expt, double p){
    return pow((1-p),no_expt)*p;
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
    //cout<<"entered" <<endl;
    vector<Node*> result=z_final;
    vec points=this->getPointIndices(z_final);
    for(int i=0;i<points.size();i++){
        //if(i==76) cout<<"Updating point "<<points[i]<<endl;
        vec index(1);index.fill(points[i]);int current=static_cast<int>(points[i]);
        //if(i==76)cout<<"y coordinate "<<this->subsetY(y,index)<<endl;
        //if(i==76)cout<<"X coordinate"<<X.row(current)<<endl;
        vec alpha=this->getSampleProbsForZ(this->subsetY(y,index),X.row(current)); //checked against R
        //if(i==76)cout<<"out of get sample probs for z"<<endl;
        //cout<<"alpha:"<<alpha<<endl;
        //cout<<"X:"<<X.row(current)<<endl;
        //cout<<"E"<<z_final[i]->id<<endl;
        //cout<<"Before:"<<result[current]->id<<endl;
        result[current]=this->updateZ_onepoint_sample(alpha);
        //cout<<"Allocated to "<<result[current]->id<<endl;
    }
    return result;
}

vector<Node*> Gate::updateZ(vec y, mat X,vector<Node*> z_final,vec* q_z_helper){
    //cout<<"entered" <<endl;
    vector<Node*> result=z_final;
    vec points=this->getPointIndices(z_final);
    for(int i=0;i<points.size();i++){
        //if(i==76) cout<<"Updating point "<<points[i]<<endl;
        vec index(1);index.fill(points[i]);int current=static_cast<int>(points[i]);
        //if(i==76)cout<<"y coordinate "<<this->subsetY(y,index)<<endl;
        //if(i==76)cout<<"X coordinate"<<X.row(current)<<endl;
        vec alpha=this->getSampleProbsForZ(this->subsetY(y,index),X.row(current)); //checked against R
        //alpha.print("alphas:");
        //if(i==76)cout<<"out of get sample probs for z"<<endl;
        //if(i==76)cout<<alpha<<endl;
        //if(i==76)cout<<"Before:"<<result[current]->name<<endl;
        result[current]=this->updateZ_onepoint_sample(alpha);
        (*q_z_helper)[i]=alpha[this->whichTerminal(result[current])];;
        //if(i==76)cout<<"Allocated to "<<result[current]->name<<endl;
    }
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
        alpha[i]=log(this->getPathProb(current,X))+sum(current->expertmodel->logdensity(y,X*beta,logsigma_sq));
        }
     alpha.elem(find_nonfinite(alpha)).fill(0);
     alpha=alpha-max(alpha);
     alpha=exp(alpha)/sum(exp(alpha));
     return alpha;
 }

//Returns the probability of being assigned to each expert for all points
 mat Gate::getSampleProbsForZ_mat(vec y, mat X){ //rows are observations and columns are experts
     mat result(X.n_rows,this->countTerminals()); 
     for(int i=0;i<X.n_rows;i++){
        //cout<<"Updating point "<<i<<endl;
        vec index(1);index.fill(i);int current=static_cast<int>(i);
        mat alpha=this->getSampleProbsForZ(this->subsetY(y,index),X.row(current));
        alpha.reshape(alpha.n_cols,alpha.n_rows);    
        result.row(i)=alpha; 
     }    
       //result.print("result");
     return result;
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
            if (rnum < alpha[i]) {
                return terminals[i];
            }else{
                rnum -= alpha[i];
            }
        }
        cout<<"we should never be here update Z"<<endl;
        alpha.print("alphas:");
        //cout<<"rnum "<<rnum<<endl;
    return terminals[alpha.size()-1];
 }

//picks a node with the probabilities given in alpha

Node* Gate::pickNode(vec alpha, vector<Node*> nodes){
    if(sum(alpha)!=1){
        double sums=sum(alpha);
        alpha=alpha/sums;   
    }
    double rnum = randu();
    for(int i=0; i<alpha.size(); i++){
        if (rnum < alpha[i]) {
            return nodes[i];
        }else{
            rnum -= alpha[i];
        }
    }
        cout<<"we should never be here"<<endl;

    return nodes[alpha.size()-1];
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

 int Gate::whichTerminal(Node* node){
     vector<Node*> terminals=this->getTerminalNodes();
     vec terminalID(terminals.size());
     for(int i=0;i<terminals.size();i++) terminalID[i]=terminals[i]->idLR;
     uvec index=find(terminalID==node->idLR);
     if(index.n_rows==0){
        return -1;
     }else{
     return static_cast<int>(as_scalar(index));
     }
 }

 int Gate::whichLastGate(Node* node){
     vector<Node*> gates=this->getLastGates();
     vec gatesID(gates.size());
     for(int i=0;i<gates.size();i++) gatesID[i]=gates[i]->idLR;
     uvec index=find(gatesID==node->idLR);
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
     if(terminals.size()>1){
        for(int i=0;i<terminals.size();i++){
        vec pathprobs=this->getPathProb_mat(terminals[i],X);
        //if(i==1)cout<<dynamic_cast<Expert*>(terminals[i])->name<<endl;
        //if(i==1)cout<<dynamic_cast<Expert*>(terminals[i])->beta<<endl;
        //if(i==1)cout<<pathprobs<<endl;
        vec est=X*(dynamic_cast<Expert*>(terminals[i])->beta); //not true for GLM
        //if(i==1)cout<<"est"<<est<<endl;
        helper.col(i)=pathprobs%est;
        }
        //helper.print("helper:");
        vec final=sum(helper,1);
        //final.print("summed up:");
        return sum(helper,1);
     }else{
         return X*(dynamic_cast<Expert*>(terminals[0])->beta);
     }
 }

vec Gate::predictions_var(mat X){
    vec preds=this->predict(X);
    vec result(preds.size()); result.fill(-1111);
    vector<Node*> terminals=this->getTerminalNodes();
    for(int i=0;i<preds.size();i++){
        vec choose(terminals.size());
        for(int j=0; j<terminals.size(); j++){
            choose[j]=this->getPathProb(terminals[j],X.row(i));
            //cout<<"Probability of choosing E"<<terminals[j]->id<<" is "<<choose[j]<<endl;
        }
        Expert* chosen=dynamic_cast<Expert*>(this->updateZ_onepoint_sample(choose));
        //cout<<"E"<<chosen->id<<" has been chosen."<<endl;
        //cout<<"Predicted value "<<preds[i]<<endl;
        //cout<<"Sigma sqared "<<exp(chosen->logsigma_sq)<<endl;
        //double helper=as_scalar(preds[i]+randn(1)*sqrt(exp(chosen->logsigma_sq))); //sqrt 
        double helper=as_scalar(preds[i]+randn(1)*exp(chosen->logsigma_sq)); //no sqrt
        result[i]=helper;
        //cout<<"Prediction with noise: "<<result[i]<<endl;
    }
    return result;
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
    //cout<<"z update called from "<<this->name<<"with ID "<<this->id<<endl;
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
        current->beta=current->expertmodel->findBeta(current->subsetY(y,current->getPointIndices(z_assign)),current->subsetX(X,current->getPointIndices(z_assign)),exp(logsigma_sq),mu_beta,Sigma_beta);
        //current->beta=current->expertmodel->findBetaMLE(current->subsetY(y,current->getPointIndices(z_assign)),current->subsetX(X,current->getPointIndices(z_assign)));
        
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
    if(this->countGates()>1){
        vector<Node*> desc=this->getDescendants();
        for(int i=0;i<desc.size();i++){
            if(desc[i]->countChildren()!=0){
                Gate* current=dynamic_cast<Gate*>(desc[i]);
                mat z=current->getZ(z_assign); 
                current->gamma=current->findGammaMLE(current->subsetX(X,current->getPointIndices(z_assign)),z,Omega);
            }
        }
    }
}

double Gate::dnorm(double y, double mu, double sigma_sq){
     return exp(-0.5*(log(2*M_PI)+log(sigma_sq))-pow(y-mu,2)/(2*sigma_sq));
}
    
void Gate::splitEmptyExpert(Expert* ExpertToSplit, Expert* ExpertToAdd, Gate* GateToAdd, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega, double* prior_betastar, double* prior_sigmastar, double* prior_gamma, double* q_betastar,double* q_sigmastar, double* q_gamma, string penalty){
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
    //0) Check the dimension of X
    int p=X_sub.n_cols-1;
    //cout<<"Dimension="<<p<<endl;
    //0.1) Project onto 1D
    mat X_transformed=X_sub;
    if(p>1){
        //cout<<"Multidimensional case"<<endl;
        vec direction(p,fill::randn);
        //cout<<"Chosen direction vector:"<<direction<<endl;
        double magnitude=static_cast<double>(as_scalar(sqrt(sum(pow(direction,2)))));
        //cout<<"Magnitude of the vector:"<<magnitude<<endl;
        direction=direction/magnitude;
        //cout<<"Standardised direction vector:"<<direction<<endl;
        mat direction_mat(p,1); direction_mat.col(0)=direction;
        //direction_mat.print("dir mat:");
        //cout<<"X_sub.n_rows"<<X_sub.n_rows<<endl;
        mat X_helper=X_sub.submat(span(0,X_sub.n_rows-1),span(1,p));
        //cout<<"X_helper: "<<X_helper.n_rows<<"x"<<X_helper.n_cols<<endl;
        mat transformed=X_helper*direction;
        //cout<<"Transformed: "<<transformed.n_rows<<"x"<<transformed.n_cols<<endl;
        mat intercept(X_sub.n_rows,1); intercept.fill(1);
        X_transformed=join_rows(intercept,transformed);
        //cout<<X_transformed.n_rows<<"x"<<X_transformed.n_cols<<endl;
    }
    //1) Draw one point at random
    int n_rand=rand() % y_sub.size();  
    mat a_star=X_transformed.row(n_rand); a_star.shed_col(0);
    //cout<<"Selected a*: "<<a_star<<endl;
    mat x_star=X_sub.row(n_rand); x_star.shed_col(0);
    //cout<<"Selected x*: "<<a_star<<endl;
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
    //gamma.print("Gamma for new gate:");
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
        //cout<<z_assign[current]->id<<endl;
        z_assign[current]=GateToAdd->updateZ_onepoint_sample(alpha);
        //cout<<z_assign[current]->id<<endl;
        (*q_z_helper)[i]=alpha[GateToAdd->whichChild(z_assign[current])];
    }

    return z_assign;
}

void Gate::proposeNewExpertParams(vec y, mat X, Expert* myExpert,vec mu_beta, mat Sigma_beta, double a, double b, double* q_betastar, double* q_sigmastar){
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    if(X.n_rows>2){
        //1) If not empty     
        //2) Calculate estimates for beta and sigma
        vec betahat;
        mat Sigma_prop;
        double sigmahat;
           if((X.t()*X).is_sympd()==0){ //if not positive definite 
            //cout<<"Not positive definite"<<endl;
            vec diags(X.n_cols); diags.fill(0.0001); diags[0]=0;
            mat varmat=diagmat(diags);
            mat X_posdef=X+(mvnrnd(mu_zeros,varmat,X.n_rows)).t();
            betahat=myExpert->expertmodel->findBetaMLE(y,X_posdef);
            sigmahat=myExpert->expertmodel->findLogSigmaSqMLE(y,X_posdef,betahat);
            Sigma_prop= (X_posdef.t()*X_posdef).i()*exp(sigmahat);
        }else{
            betahat=myExpert->expertmodel->findBetaMLE(y,X);
            sigmahat=myExpert->expertmodel->findLogSigmaSqMLE(y,X,betahat);
            Sigma_prop=(X.t()*X).i()*exp(sigmahat);
        } 
        //betahat.print("betahat:");
        //cout<<"sigmahat:"<<sigmahat<<endl;
        //Sigma_prop.print("Sigma_prop:");
        //3) Draw a value for beta centred around betahat
        myExpert->beta=betahat+mvnrnd(mu_zeros, Sigma_prop,1);
        //cout<<"Proposed beta: "<<myExpert->beta<<endl;
        //4) Draw a value for sigma given beta
        myExpert->logsigma_sq=myExpert->expertmodel->updateSigma(y,X,myExpert->beta,a,b,X.n_rows);
        //cout<<"Proposed sigma: "<< myExpert->logsigma_sq<<endl;
        //5) Record the densities
        *q_betastar+=sum(this->logmvndensity(myExpert->beta,betahat,Sigma_prop));
        *q_sigmastar+=myExpert->expertmodel->qSigma(y,X,myExpert->beta,myExpert->logsigma_sq,a,b);
        //cout<<"q_sigma: "<<*q_sigmastar<<endl;
        //cout<<"a+X.n_rows/2: "<<a+static_cast<double>(X.n_rows)/2<<endl;
    }else{
        //1) If empty
        //2) Draw both beta and sigma from prior
        myExpert->beta=mvnrnd(mu_beta, Sigma_beta,1);
        myExpert->logsigma_sq=log(1/randg( distr_param(a,1/b)));
        *q_betastar+=sum(this->logmvndensity(myExpert->beta,mu_beta,Sigma_beta));
        *q_sigmastar+=myExpert->expertmodel->qSigma(y,X,myExpert->beta,myExpert->logsigma_sq,a,b);
    }
    //cout<<"beta:"<<myExpert->beta<<endl;
    //cout<<"sigma:"<<myExpert->logsigma_sq<<endl;
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

vector<Node*> Gate::split(vec y, mat X, int* accept, Expert* ExpertToSplit, Expert* ExpertToAdd, Gate* GateToAdd, vector<Node*> z_assign,  vec mu_beta, mat Sigma_beta, vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, double a, double b, mat Omega, double lambda, string penalty){// double mu_jump, double sigma_jump, vec mu_beta, mat Sigma_beta, vec* x_record, vec* gamma_record){
    //1) Record the indices of points that reached the expert we are trying to split
    vec  points=ExpertToSplit->getPointIndices(z_assign);
    //cout<<"Number of points in the expert to split: "<<points.size()<<endl;
    //2) Check if the chosen expert is empty
    bool empty=points.size()<=2;
    //cout<<"Empty: "<<empty<<endl;
    //3) Subset y and X accordingly
    vec y_sub=ExpertToSplit->subsetY(y,points);
    mat X_sub=ExpertToSplit->subsetX(X,points);
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    //4) Record relevant information before splitting
    double q_beta=ExpertToSplit->expertmodel->qBeta(y_sub,X_sub,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,mu_beta,Sigma_beta);
    double q_sigma=ExpertToSplit->expertmodel->qSigma(y_sub,X_sub,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,a,b);
    double prior_beta=sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    double prior_sigma=ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    double size_prior=0;
    if(penalty=="poi") size_prior=log(this->poissondensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
    if(penalty=="geo") size_prior=log(this->geomdensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
    double p_Model=log(ExpertToSplit->mostSeniorGate()->prob_RJ(z_assign,"split")[ExpertToSplit->mostSeniorGate()->whichTerminal(ExpertToSplit)]);
    //cout<<"E"<<ExpertToSplit->id<<" was chosen to split with prob "<<exp(p_Model)<<endl;
    //5) If chosen expert isn't empty, record its log-likelihood
    double loglik_old=100000000;
    if(empty==0) loglik_old=ExpertToSplit->expertmodel->loglik(y_sub,X_sub*(ExpertToSplit->beta),ExpertToSplit->logsigma_sq);
    
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
    double size_prior_new=0;
    vector<Node*> z_new(z_assign.size());

    //9)Make parameter proposals
    if(empty==1){
        cout<<"Splitting an empty expert"<<endl;
        this->splitEmptyExpert(ExpertToSplit,ExpertToAdd, GateToAdd, mu_beta,Sigma_beta,a,b,Omega,&prior_betastar,&prior_sigmastar,&prior_gamma,&q_betastar,&q_sigmastar,&q_gamma,penalty);
        double size_prior_new=0;
        z_new=z_backup;
        double p_Model_new=log(GateToAdd->mostSeniorGate()->prob_RJ(z_new,"merge")[ExpertToSplit->mostSeniorGate()->whichLastGate(GateToAdd)]);  
        if(penalty=="poi") size_prior_new = log(this->poissondensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
        if(penalty=="geo") size_prior_new = log(this->geomdensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
        acceptance=size_prior_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma+p_Model_new-
                   size_prior-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_gamma-p_Model;
    }else{
        //cout<<"Splitting a full expert"<<endl;
        //10) Propose a gamma
        GateToAdd->gamma=GateToAdd->proposeGammaSplit(y_sub,X_sub,mu_gamma1,Sigma_gamma1,sigma_epsilon);
        //11) Propose allocations after the split
        vec q_z_helper(y_sub.size());
        z_new=this->proposeZafterSplit(y_sub,X_sub, points, GateToAdd,z_assign,&q_z_helper);
        //12) Subset accordingly
        vec y1=ExpertToAdd->subsetY(y,ExpertToAdd->getPointIndices(z_new)); mat X1=ExpertToAdd->subsetX(X,ExpertToAdd->getPointIndices(z_new)); int n1=static_cast<int>(X1.n_rows);
        vec y2=ExpertToSplit->subsetY(y,ExpertToSplit->getPointIndices(z_new)); mat X2=ExpertToSplit->subsetX(X,ExpertToSplit->getPointIndices(z_new));int n2=static_cast<int>(X2.n_rows);
        //13) Propose expert parameters after the split
        this->proposeNewExpertParams(y1,X1,ExpertToAdd,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
        this->proposeNewExpertParams(y2,X2,ExpertToSplit,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
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
        //18) Calculate model size prior
        if(penalty=="poi") size_prior_new = log(this->poissondensity(GateToAdd->mostSeniorGate()->countTerminals(),lambda));
        if(penalty=="geo") size_prior_new = log(this->geomdensity(GateToAdd->mostSeniorGate()->countTerminals(),lambda));
        //19) Calculate the probability of going back from the split expert to it being merged
        double p_Model_new=log(GateToAdd->mostSeniorGate()->prob_RJ(z_new,"merge")[ExpertToSplit->mostSeniorGate()->whichLastGate(GateToAdd)]);
        //cout<<"G"<<GateToAdd->id<<" would be selected to merge with prob "<<exp(p_Model_new)<<endl;
        //20) Calculate acceptance
        // //acceptance=loglik_new+size_prior_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma+p_Model_new-
        //   //         loglik_old-size_prior-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_z-q_gamma-p_Model;//a-q_z2; 
        // cout<<"q_z"<<q_z<<endl;
        acceptance=loglik_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma+p_Model_new-
                   loglik_old-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_z-q_gamma-p_Model;
    }


   // UNCOMMENT BELOW IF WANT TO ZOOM IN ON VALUES
    // cout<<"SPLIT"<<endl;
    // cout<<"acceptance"<<acceptance<<endl;
    // cout<<"loglik_old"<<loglik_old<<endl;
    // cout<<"loglik_new"<<loglik_new<<endl;
    // cout<<"prior_betastar"<<prior_betastar<<endl;
    // cout<<"prior_beta"<<prior_beta<<endl;
    // cout<<"prior_sigmastar"<<prior_sigmastar<<endl;
    // cout<<"prior_sigma"<<prior_sigma<<endl;
    // cout<<"prior_gamma"<<prior_gamma<<endl;
    // cout<<"q_betastar"<<q_betastar<<endl;
    // cout<<"q_beta"<<q_beta<<endl;
    // cout<<"q_sigmastar"<<q_sigmastar<<endl;
    // cout<<"q_sigma"<<q_sigma<<endl;
    // cout<<"q_gamma:"<<GateToAdd->gamma<<endl;
    // cout<<"q_gamma"<<q_gamma<<endl;
    // cout<<"size_prior"<<size_prior<<endl;
    // cout<<"size_prior_new"<<size_prior_new<<endl;


    double u=randu();
    bool acc=log(u)<=acceptance;
    if(acc==1){
        //cout<<"Split has been accepted"<<endl;
        // cout<<ExpertToSplit->name<<" split into "<< GateToAdd->Children[0]->name<<" and "<< GateToAdd->Children[1]->name<<endl;
        // vector<Node*> gates=GateToAdd->mostSeniorGate()->getGates();
        // vector<Node*> experts=GateToAdd->mostSeniorGate()->getTerminalNodes();
        // for(int i=0;i<gates.size();i++) dynamic_cast<Gate*>(gates[i])->printChildren();
        // for(int j=0;j<experts.size();j++) cout<<experts[j]->name<<" has id "<<experts[j]->id<<endl;
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


vector<Node*> Gate::split_at_root(vec y, mat X, int* accept, Expert* ExpertToSplit, Expert* ExpertToAdd, vector<Node*> z_assign,  vec mu_beta, mat Sigma_beta, vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, double a, double b, mat Omega, double lambda, string penalty){// double mu_jump, double sigma_jump, vec mu_beta, mat Sigma_beta, vec* x_record, vec* gamma_record){
    //1) Record the indices of all points as integers
    vec  points=ExpertToSplit->getPointIndices(z_assign);
    //2) Set up helpers
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    //3) Record relevant information before splitting
    double q_beta=ExpertToSplit->expertmodel->qBeta(y,X,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,mu_beta,Sigma_beta);
    double q_sigma=ExpertToSplit->expertmodel->qSigma(y,X,ExpertToSplit->beta,ExpertToSplit->logsigma_sq,a,b);
    double prior_beta=sum(this->logmvndensity(ExpertToSplit->beta,mu_beta,Sigma_beta));
    double prior_sigma=ExpertToSplit->expertmodel->IG_log(exp(ExpertToSplit->logsigma_sq),a,b);
    double size_prior=0;
    if(penalty=="poi") size_prior=log(this->poissondensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
    if(penalty=="geo") size_prior=log(this->geomdensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
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
    this->proposeNewExpertParams(y1,X1,ExpertToAdd,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    this->proposeNewExpertParams(y2,X2,ExpertToSplit,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
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
    //16) Calculate size prior
    double size_prior_new=0;
    if(penalty=="poi") size_prior_new=log(this->poissondensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
    if(penalty=="geo") size_prior_new=log(this->geomdensity(ExpertToSplit->mostSeniorGate()->countTerminals(),lambda));
    //16) Calculate acceptance
    acceptance=loglik_new+size_prior_new+prior_betastar+prior_sigmastar+prior_gamma+q_beta+q_sigma-
               loglik_old-size_prior-prior_beta-prior_sigma-q_betastar-q_sigmastar-q_z-q_gamma; 
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


vector<Node*> Gate::merge(vec y, mat X, int* accept, Gate* GateToMerge, vector<Node*> z_assign, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega, double lambda, string penalty){
    vec mu_zeros(X.n_cols); mu_zeros.fill(0);
    //1) Determine which experts to merge
    vector<Node*> ExpertsToMerge=GateToMerge->getTerminalNodes();
    //cout<<"E"<<ExpertsToMerge[0]->id<<" and E"<<ExpertsToMerge[1]->id<<" will be merged."<<endl;
    //2) Subset the points that have reached the gate
    vec points=GateToMerge->getPointIndices(z_assign);
    //points.print("points");
    //cout<<"Total points: "<<points.size()<<endl;
    //3) Record the log likelihood, priors and parameter densities for the current architecture
    double loglik_old=GateToMerge->loglik_complete(y,X,z_assign); 
    double prior_beta=0;
    double prior_sigma=0;
    double q_beta=0;
    double q_sigma=0;
    double q_gamma=GateToMerge->q_gammaMerge(X,z_assign,Omega); //subsets X based on z_assign inside and deals with an empty case
    double prior_gamma=sum(this->logmvndensity(GateToMerge->gamma,mu_zeros,Omega)); 
    double size_prior=0;
    double p_Model=log(GateToMerge->mostSeniorGate()->prob_RJ(z_assign,"merge")[GateToMerge->mostSeniorGate()->whichLastGate(GateToMerge)]);
    //cout<<"G"<<GateToMerge->id<<" was chosen to be merged with probability of "<<exp(p_Model)<<endl;
    if(penalty=="poi") size_prior=log(this->poissondensity(GateToMerge->mostSeniorGate()->countTerminals(),lambda));
    if(penalty=="geo") size_prior=log(this->geomdensity(GateToMerge->mostSeniorGate()->countTerminals(),lambda));
    for(int i=0;i<ExpertsToMerge.size();i++){
        Expert* current=dynamic_cast<Expert*>(ExpertsToMerge[i]);
        vec points_helper=current->getPointIndices(z_assign);
        vec y_sub=current->subsetY(y,points_helper);
        mat X_sub=current->subsetX(X,points_helper);
        prior_beta=prior_beta+sum(this->logmvndensity(current->beta,mu_beta,Sigma_beta));
        prior_sigma=prior_sigma+current->expertmodel->IG_log(exp(current->logsigma_sq),a,b);
        //cout<<"Sigma:"<<exp(current->logsigma_sq)<<endl;
        //cout<<"Prior sigma for E"<<current->id<<": "<<prior_sigma<<endl;
        q_beta=q_beta+current->expertmodel->qBeta(y_sub,X_sub,current->beta,current->logsigma_sq,mu_beta,Sigma_beta);
        q_sigma=q_sigma+current->expertmodel->qSigma(y_sub,X_sub,current->beta,current->logsigma_sq,a,b);
        //cout<<"q_sigma merge: "<<q_sigma<<endl;
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
        //cout<<"Merging gate G"<<GateToMerge->id<<endl;
        GateToMerge->Parent->replaceChild(k,ExpertsToMerge[0]);
        //cout<<"The parent of E"<<ExpertsToMerge[0]->id<<" now is G"<<ExpertsToMerge[0]->Parent->id<<endl;
        //cout<<"Another child of G"<<ExpertsToMerge[0]->Parent->id<<" is G"<<ExpertsToMerge[0]->Parent->Children[0]->id<<endl;
        //for(int i=0; i<ExpertsToMerge[0]->Parent->Children[0]->countChildren();i++) cout<<"G"<<ExpertsToMerge[0]->Parent->Children[0]->id<<" has a child E"<<ExpertsToMerge[0]->Parent->Children[0]->Children[i]->id<<endl;
        
    }
    //7) Create a new pointer to the newly formed expert to simplify notation
    Expert* MergedExpert=dynamic_cast<Expert*>(ExpertsToMerge[0]);
    //MergedExpert->mostSeniorGate()->printChildren();
    //8) Issue IDs to the newly created architecture
    MergedExpert->mostSeniorGate()->issueID();
    MergedExpert->mostSeniorGate()->issueIDLR();
    vector<Node*> terminals=MergedExpert->mostSeniorGate()->getTerminalNodes();
    //for(int i=0;i<terminals.size();i++) cout<<"E"<<terminals[i]->id<< " parent is G"<<terminals[i]->Parent->id<<endl;
    
    //9) Assign all points to one merged expert
    for(int i=0;i<points.size();i++){
         z_new[static_cast<int>(points[i])]=MergedExpert;
    }
    //10) Subset all points that are in the new merged expert
    mat myX=MergedExpert->subsetX(X,MergedExpert->getPointIndices(z_new));
    vec myY=MergedExpert->subsetY(y,MergedExpert->getPointIndices(z_new));
    //11) Propose post merge parameters
    double q_betastar=0; double q_sigmastar=0;
    MergedExpert->Parent->proposeNewExpertParams(myY,myX,MergedExpert,mu_beta,Sigma_beta,a,b,&q_betastar,&q_sigmastar);
    //cout<<"Merged Expert Params:"<<endl;
    //cout<<"beta: "<<MergedExpert->beta<<endl;
    //cout<<"sigma: "<<MergedExpert->logsigma_sq<<endl;
    //11*) Perform one update of the allocations
    vec q_z_helper(y.size());
    z_new=MergedExpert->mostSeniorGate()->updateZ(y,X,z_new,&q_z_helper);
    double q_z=sum(log(q_z_helper));
    //12) Record log-likelihood after the merge
    double loglik_new=MergedExpert->expertmodel->loglik(myY,myX*(MergedExpert->beta),MergedExpert->logsigma_sq);
    //13) Evaluate priors
    double prior_betastar=sum(this->logmvndensity(MergedExpert->beta,mu_beta,Sigma_beta));
    double prior_sigmastar=MergedExpert->expertmodel->IG_log(exp(MergedExpert->logsigma_sq),a,b);
    //cout<<"sigma star: "<<exp(MergedExpert->logsigma_sq)<<endl;
    //cout<<"prior sigmastar: "<<prior_sigmastar<<endl;
    //14) Calculate size prior
    double size_prior_new=0;
    if(penalty=="poi") size_prior_new=log(this->poissondensity(MergedExpert->mostSeniorGate()->countTerminals(),lambda));
    if(penalty=="geo") size_prior_new=log(this->geomdensity(MergedExpert->mostSeniorGate()->countTerminals(),lambda));
    //15) Calculate the probability of choosing to split the merged expert
    double p_Model_new=log(MergedExpert->mostSeniorGate()->prob_RJ(z_new,"split")[MergedExpert->mostSeniorGate()->whichTerminal(MergedExpert)]);
    //cout<<"E"<<MergedExpert->id<<" would be chosen to split with probability of "<<exp(p_Model_new)<<endl;
    //15) Calculate the acceptance probability
    //double acceptance=loglik_new+size_prior_new+prior_betastar+prior_sigmastar+q_beta+q_sigma+q_gamma+p_Model_new-
      //               loglik_old-size_prior-prior_beta-prior_sigma-prior_gamma-q_betastar-q_sigmastar-q_z-p_Model;
     double acceptance=loglik_new+prior_betastar+prior_sigmastar+q_beta+q_sigma+q_gamma+p_Model_new-
                       loglik_old-prior_beta-prior_sigma-prior_gamma-q_betastar-q_sigmastar-q_z-p_Model;                  
  

    //UNCOMMENT IF WANT TO SEE
    // cout<<"MERGE"<<endl;
    // cout<<"acceptance"<<acceptance<<endl;
    // cout<<"loglik_old"<<loglik_old<<endl;
    // cout<<"loglik_new"<<loglik_new<<endl;
    // cout<<"prior_betastar"<<prior_betastar<<endl;
    // cout<<"prior_beta"<<prior_beta<<endl;
    // cout<<"prior_sigmastar"<<prior_sigmastar<<endl;
    // cout<<"prior_sigma"<<prior_sigma<<endl;
    // cout<<"prior_gamma"<<prior_gamma<<endl;
    // cout<<"q_betastar"<<q_betastar<<endl;
    // cout<<"q_beta"<<q_beta<<endl;
    // cout<<"q_sigmastar"<<q_sigmastar<<endl;
    // cout<<"q_sigma"<<q_sigma<<endl;
    // //cout<<"gamma:"<<GateToMerge->gamma<<endl;
    // cout<<"q_gamma"<<q_gamma<<endl;
    // //cout<<"size_prior"<<size_prior<<endl;
    // //cout<<"size_prior_new"<<size_prior_new<<endl;
    // cout<<"p_Model"<<p_Model<<endl;
    // cout<<"p_Model_new"<<p_Model_new<<endl;

   //16) Check if accepted
    double u=randu();
    acc=log(u)<acceptance;
    if(acc==1){
        //cout<<"Merge has been accepted"<<endl;
        // cout<<"The kids of "<<backup->name<<" have been merged"<<endl;
        // cout<<": "<<backup->Children[0]->name<<" and "<< backup->Children[1]->name<<endl;
        // vector<Node*> gates=MergedExpert->mostSeniorGate()->getGates();
        // vector<Node*> experts=MergedExpert->mostSeniorGate()->getTerminalNodes();
        // for(int i=0;i<gates.size();i++) dynamic_cast<Gate*>(gates[i])->printChildren();
        // for(int j=0;j<experts.size();j++) cout<<experts[j]->name<<" has id "<<experts[j]->id<<endl;
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

 double Gate::q_gammaMerge(mat X, vector<Node*> z_assign, mat Omega){
     mat z=this->getZ(z_assign);
     mat myX=this->subsetX(X,this->getPointIndices(z_assign));
     mat R;
     vec gammahat = this->findGammaQR(myX, z, Omega,&R);
     //gammahat.print("gammahat:");
     //cout<<"gamma:"<<this->gamma<<endl;
     //cout<<"MLE:"<<this->findGammaMLE(myX,z,Omega)<<endl;
     mat Sigma=(R.t()*R).i();
     //Sigma.print("Sigma:");
     //cout<<"result:"<<this->logmvndensity(this->gamma,gammahat,Sigma)<<endl;
     return sum(this->logmvndensity(this->gamma,gammahat,Sigma));
 }

vector<Node*> Gate::MCMC_RJ(int N,  bool doRJ, int RJ_every, int L, mat* accept_RJ, vec y, mat X, vec mu_beta, mat Sigma_beta, double a, double b, mat Omega,vec mu_gamma1, mat Sigma_gamma1, double sigma_epsilon, vector<Node*> z_final, mat X_new, mat* predictions, mat* predictions_var, int predict_every, int record_params_every, mat* no_expt, double lambda, mat* z_record, mat* beta_record, mat* sigma_record, mat* gamma_record, mat* pi_record, mat* pi_record2,string penalty, mat* dens){
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
            G0=dynamic_cast<Gate*>(RootParent->Children[0]); //record current root node
            if(i%100==0) cout<<"Iteration number "<<i<<endl;
            //Run MCMC, if only one expert, then only update expert parameters
            if(G0->countTerminals()>1){ 
                z_final=G0->MCMC_OneRun(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }else{
                dynamic_cast<Expert*>(G0->Children[0])->MCMC_internal(y,X,mu_beta,Sigma_beta,a,b,Omega,z_final);
            }
            //Record allocations:
            for(int h=0;h<y.size();h++){
                (*z_record).submat(span(i,i),span(h,h))=z_final[h]->id;
            }
            //Record predictions:
            if(i%predict_every==0){
                        (*predictions).col(predict_no)=G0->predict(X_new);
                        (*predictions_var).col(predict_no)=G0->predictions_var(X_new);
                        predict_no+=1;
            } //closes predict_every
            //Record betas:
            //if((i-1)%RJ_every==0){
            if(i%record_params_every==0){
            vector<Node*> terminals=G0->getTerminalNodes(); 
            mat beta_record_helper(X.n_cols+2,terminals.size()); beta_record_helper.fill(i);
            mat sigma_record_helper(3,terminals.size()); sigma_record_helper.fill(i);
            for(int g=0;g<terminals.size();g++){
                beta_record_helper.submat(span(1,X.n_cols),span(g,g))=dynamic_cast<Expert*>(terminals[g])->beta;
                beta_record_helper.submat(span(X.n_cols+1,X.n_cols+1),span(g,g))=terminals[g]->id;
                sigma_record_helper.submat(span(1,1),span(g,g))=dynamic_cast<Expert*>(terminals[g])->logsigma_sq;
                sigma_record_helper.submat(span(2,2),span(g,g))=terminals[g]->id;
            }
            *beta_record=join_rows(*beta_record,beta_record_helper);
            *sigma_record=join_rows(*sigma_record,sigma_record_helper);
            //}
            //Record gammas:
           //  if((i-1)%RJ_every==0){
            vector<Node*> gates=G0->getGates(); mat gamma_record_helper(X.n_cols+2,gates.size()); gamma_record_helper.fill(i);
            for(int h=0;h<gates.size();h++){
                gamma_record_helper.submat(span(1,X.n_cols),span(h,h))=dynamic_cast<Gate*>(gates[h])->gamma;
                gamma_record_helper.submat(span(X.n_cols+1,X.n_cols+1),span(h,h))=gates[h]->id;
            }
            *gamma_record=join_rows(*gamma_record,gamma_record_helper);
            //}
            //record pi's
           // if((i-1)%RJ_every==0){
            vector<Node*> experts=G0->getTerminalNodes(); 
            mat pi_record_helper(X.n_rows+2,experts.size()); pi_record_helper.fill(i);
            mat pi_record_helper2(X.n_rows+2,experts.size()); pi_record_helper2.fill(i);
            for(int v=0;v<experts.size();v++){
                pi_record_helper.submat(span(1,1),span(v,v))=experts[v]->id;
                pi_record_helper2.submat(span(1,1),span(v,v))=experts[v]->id;
                pi_record_helper2.submat(span(2,2+X.n_rows-1),span(v,v))=G0->getPathProb_mat(experts[v],X);
            }
            pi_record_helper.submat(span(2,2+X.n_rows-1),span(0,experts.size()-1))=G0->getSampleProbsForZ_mat(y,X);
            *pi_record=join_rows(*pi_record,pi_record_helper);
            *pi_record2=join_rows(*pi_record2,pi_record_helper2);
            }
            
            if(doRJ==1){ //If the RJ is on
                    if(i%RJ_every==0){ //Check if it is the correct iteration
                    //cout<<"Jump number "<<jump_no<<endl;
                        for(int j=0;j<L;j++){   
                            G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                            //cout<<"RJ taking place at iteration "<< i<<" for the "<<j+1<<"-th time"<<endl;
                            int RJ_direction= rand() % 2; //Randomly choose the direction of the jump
                            if(G0->mostSeniorGate()->countTerminals()==1) RJ_direction=0; //If there is only one expert propose to split it
                            if(RJ_direction==0){ //in case of a split
                                //cout<<"Performing a split."<<endl;                            
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);                               
                                vector<Node*> terminals=G0->getTerminalNodes();                                
                                vec p_split=G0->prob_RJ(z_final,"split");                               
                                Expert* ExpertToSplit = dynamic_cast<Expert*>(this->pickNode(p_split,terminals));
                                //cout<<"Chosen to split E"<<ExpertToSplit->id<<endl;
                                //Commented out proposing to split expert at random
                                //int n_rand=rand() % terminals.size(); //chose expert to split at random
                                //Expert* ExpertToSplit = dynamic_cast<Expert*>(terminals[n_rand]);

                                Expert* ExpertToAdd= new Expert(); //create a new sibling expert
                                ExpertToAdd->name= "E"+to_string(G0->getMaxExpertID()+1);
                                NormalFamily* NF=new NormalFamily();
                                ExpertToAdd->expertmodel=NF;
                                int accept;
                                if(ExpertToSplit->Parent->countChildren()==1){
                                    //cout<<"Splitting at Root"<<endl;
                                    z_final=G0->split_at_root(y,X,&accept,ExpertToSplit,ExpertToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega,lambda,penalty);
                                }else{
                                    Gate*   GateToAdd=new Gate(); //create a new gate parent
                                    GateToAdd->name="G"+to_string(G0->getMaxGateID()+1);
                                    //G0->printChildren();
                                    z_final=G0->split(y,X,&accept,ExpertToSplit,ExpertToAdd,GateToAdd,z_final,mu_beta,Sigma_beta,mu_gamma1,Sigma_gamma1,sigma_epsilon,a,b,Omega,lambda,penalty);
                                }//closes the usual split when more than 1 child
                                //RootParent->Children[0]->printChildren();
                                //ExpertToSplit->Parent->printChildren();
                                (*accept_RJ).submat(span(jump_no,jump_no),span(0,0))=1;
                                (*accept_RJ).submat(span(jump_no,jump_no),span(1,1))=accept;
                                (*no_expt).submat(span(jump_no,jump_no),span(0,0))=dynamic_cast<Gate*>(RootParent->Children[0])->countTerminals();
                                (*no_expt).submat(span(jump_no,jump_no),span(1,1))=dynamic_cast<Gate*>(RootParent->Children[0])->countEmptyTerminals(z_final);
                                jump_no+=1;
                                //cout<<"Post split"<<endl;
                                //vector<Node*> gates0=dynamic_cast<Gate*>(RootParent->Children[0])->getGates();
                                //for(int m=0;m<gates0.size();m++) gates0[m]->printChildren();
                            }else{
                                //cout<<"Performing a merge"<<endl;
                                G0=dynamic_cast<Gate*>(RootParent->Children[0]);
                                vector<Node*> gates=G0->getLastGates();
                                //for(int m=0;m<gates0.size();m++) gates0[m]->printChildren();//cout<<gates[m]->name<<endl;
                                vec p_merge=G0->prob_RJ(z_final,"merge");
                                Gate* GateToMerge;
                                GateToMerge=dynamic_cast<Gate*>(this->pickNode(p_merge,gates));                                                         
                                //Commented out choosing gate to merge at random:
                                //int n_rand=rand() % gates.size();
                                //GateToMerge=dynamic_cast<Gate*>(gates[n_rand]); //choose a gate to merge at random
                                //cout<<"Propose to merge the children of G"<<GateToMerge->id<<endl;
                                //GateToMerge->printChildren();
                                //vector<Node*> last=G0->getTerminalNodes();
                                //for(int v=0;v<last.size();v++) cout<<last[v]->name<<" has id "<<last[v]->id<<endl;
                                //for(int h=0;h<GateToMerge->countChildren();h++) cout<<GateToMerge->Children[h]->name<<" has ID "<<GateToMerge->Children[h]->id<<endl;
                                int accept2;
                                z_final=G0->merge(y,X,&accept2,GateToMerge,z_final,mu_beta,Sigma_beta,a,b,Omega,lambda,penalty);
                                //RootParent->Children[0]->printChildren();
                                (*accept_RJ).submat(span(jump_no,jump_no),span(0,0))=2;
                                (*accept_RJ).submat(span(jump_no,jump_no),span(1,1))=accept2;
                                (*no_expt).submat(span(jump_no,jump_no),span(0,0))=dynamic_cast<Gate*>(RootParent->Children[0])->countTerminals();
                                (*no_expt).submat(span(jump_no,jump_no),span(1,1))=dynamic_cast<Gate*>(RootParent->Children[0])->countEmptyTerminals(z_final);
                                jump_no+=1;
                                //cout<<"Post Merge"<<endl;
                                //vector<Node*> gates2=dynamic_cast<Gate*>(RootParent->Children[0])->getGates();
                                //for(int m=0;m<gates2.size();m++) gates2[m]->printChildren();//cout<<gates[m]->name<<endl;
                                
                            } //closes the merge
                        } //closes for L
                    } //closes if divisible    
                } //closes if doRJ
                } //closes the first for loop
            return z_final;
        } //closes the function

 vec Gate::prob_RJ(vector<Node*> z_final, string direction){
     vec result;
     if(direction=="split"){
         vector<Node*> terminals = this->getTerminalNodes();
         result.set_size(terminals.size());
         for(int p=0;p<terminals.size();p++){
             //cout<<"Calculation for a split"<<endl;
              result[p]=static_cast<double>(terminals[p]->countPoints(z_final)+1e-6)/z_final.size();
              //cout<<"E"<<terminals[p]->id<<" "<<result[p]<<endl;
         }
     }else{
         //cout<<"Calculation for a merge"<<endl;
         vector<Node*> gates=this->getLastGates();
         result.set_size(gates.size());
         for(int p=0; p<gates.size(); p++){
              result[p]=1/static_cast<double>(gates[p]->countPoints(z_final)+1e-6);
              //cout<<"G"<<gates[p]->id<<" "<<result[p]<<endl;
         }
    }
    result=result/sum(result);
    //result.print("standardised:");
    return result;
 }
   

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