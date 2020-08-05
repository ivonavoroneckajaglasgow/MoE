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

Gate::Gate(){
    cout<<"Gate has been created."<<endl;
}

void Gate::addChild(Node* aChild){
    this -> Children.push_back(aChild);
    aChild -> Parent = this;
    cout<<"Child has been added to the parent."<<endl;
}

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

void Gate::printDescendants(){

    for(int i=0;i<Children.size();i++){
        cout<<Children[i]->name<<endl;
        if(Children[i]->Children.size()!=0){
            Children[i]->printDescendants();
        }
    }
}

void Gate::printTerminalNodes(){
    
    for (int i = 0; i < Children.size(); i++) {
        if (Children[i]->Children.size() == 0) {
            cout << Children[i]->name << endl;
        } else {
            Children[i]->printTerminalNodes();
        }
    }
}
vector<Node*> Gate::getChildren() {
    return Children;
}

/**
 * @brief An internal function which helps retrieve all descendents of the gate.
 * An internal function that is then called at the Node level to retrieve all descendants of the gate.
 * @param desc vector to be filled in with the descendants of the gate.
 * @return vector<Node*> vector of pointers to the descendants of the gate.
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
 * @brief An internal function which helps retrieve all terminal nodes descending from the gate.
 * An internal function that is then called at the node level to retrieve all terminal nodes of the gate.
 * @param terminal vector to be filled in with the terminal nodes descending from the gate.
 * @return vector<Node*> vector of pointers to the terminal nodes descending from the gate.
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
};

int Gate::countChildren(){
    int n;
    n=Children.size();
    return n;
};

/**
 * @brief Counts the number of descendants of the gate.
 * 
 * @return int integer number of descendants of the gate.
 */
int Gate::countDescendants(){
    vector<Node*> desc;
    desc=this->getDescendants();
    return desc.size();
};

int Gate::refId(){
    return Children[0]->id;
}

/**
 * @brief A top layer function for assigning ID.
 * 
 * Sets the gate id to start at 2 (the root gate is automatically determined in issueID_helper2()) and 
 * the expert id to start at 1. Calls the helper functions to perform the task.
 * 
 */
void Gate::issueID(){

    int gateid=2;
    int expertid=1;

    this->issueID_helper2(&gateid,&expertid);
}
/**
 * @brief First internal function which helps to assign the id's to the nodes of tree vertically.
 * This function checks each child of the Node. If a child doesn't have any children, it assumes that it is an expert and 
 * assigns expert id to it. If a child has some children, a gate id is assigned. 
 * @param gate_id pointer to an integer which tracks gate id's.
 * @param expert_id pointer to an integer which tracks expert id's.
 */
void Gate::issueID_helper1(int* gate_id, int* expert_id){

    for(int i=0;i<this->Children.size();i++){
        if(this->Children[i]->countChildren()==0){
            cout<<"I know "<< this->Children[i]->name <<" is an expert so I assign the ID "<<*expert_id<<endl;
            this->Children[i]->id=*expert_id;
            (*expert_id)++;
        }else{
            cout<<"I know "<< this->Children[i]->name <<" is a gate so I assign ID "<<*gate_id<<endl;
            this->Children[i]->id=*gate_id;
            (*gate_id)++;
        }
    }
}

/**
 * @brief Second internal function which helps to assign the id's to the nodes of tree vertically.
 * This function identifies a root gate by checking the presence of the parent node. If there is no parent 
 * node, an id of 1 is assigned. Next, issueID_helper1() is called to assign the id's to the children of the 
 * gate. The process is repeated for every node in the tree by calling this function recursively.
 * @param gate_id pointer to an integer which tracks gate id's.
 * @param expert_id pointer to an integer which tracks expert id's.
 */
void Gate::issueID_helper2(int* gate_id, int* expert_id){

    if(this->Parent==NULL){
        this->id=1;
        cout<<"I know "<<this->name<<" is a root gate, so I assign ID "<<this->id<<endl;
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

double Gate::loglik(mat z, mat pi){
    return sum(vectorise(z%log(pi)))+sum((1-sum(z,1))%log(1-sum(pi,1)));
}

double Gate::loglik(mat X, vec gamma, mat z){
    int p=X.n_cols;
    int r=z.n_cols;
    mat gamma2=gamma;
    gamma2.reshape(p,r);
    return sum(sum(z%(X*gamma2),1)-log(1+sum(exp(X*gamma2),1)));
}

mat Gate::pi_calculator(mat X, vec gamma){
    int p=X.n_cols;
    int r=gamma.size()/p;
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
    //return arma::min(1-E,arma::max(final,E));
}

mat Gate::getRowSumsMat(mat A){
    mat rowsums(A.n_rows,A.n_cols);
    for(int i=0; i<A.n_cols; i++){
        rowsums.col(i)=sum(A,1);
    }
    return rowsums;
}

 mat Gate::pi_calculator2(mat X, vec gamma) {
    int p = X.n_cols;
    int rp= gamma.size();
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

mat Gate::pi_calculator3(mat X, vec gamma) {
    int p = X.n_cols;
    int rp= gamma.size();
    int r = rp/p;
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

vec Gate::score(mat X, mat z, mat pi){
    return vectorise(X.t()*(z-pi));
}

mat Gate::hessian(mat X, mat pi){
 int r=pi.n_cols;
 int p=X.n_cols;

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

vec Gate::findGammaMLEChol(mat X, mat z, mat Omega){
    mat L;
    return this->findGammaMLEChol(X,z,Omega,&L);
}

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

vec Gate::findGammaMLEQR(mat X, mat z, mat Omega){
    mat R;
    return this->findGammaMLEQR(X,z,Omega,&R);
}


vec Gate::findGammaMLE(mat X, mat z, mat Omega){
    vec gamma(X.n_cols*z.n_cols);
    gamma.zeros();
    for(int i=0; i<100; i++){
        mat pi=this->pi_calculator(X,gamma);
        vec gamma_old=gamma;
        gamma=gamma-solve(this->hessian(X,pi)-Omega,this->score(X,z,pi)-Omega*gamma);
    if(all(abs(gamma-gamma_old)<(sqrt(EPS),sqrt(EPS)*abs(gamma)).max())) break;
    }
    return gamma;
}

mat Gate::makeAchol(vec pi){
  mat helper(pi.size(),pi.size());
  for(int i=0; i<pi.size();i++){
      for(int j=0; j<pi.size(); j++){
          helper(i,j)=pi[i]*pi[j];
      }
  }
  return chol(diagmat(pi)-helper);
}

mat Gate::makeAchol2(vec pi){
  mat B(pi.size(),pi.size());
  B.fill(-1);
  vec b= (1-pi)/pi;
  B.diag()=b;
  mat D=diagmat(pi);
  mat L=chol(B);
  return L*D;
}

mat Gate::getXout(mat X, mat pi){
    int r=pi.n_cols;
    int n=X.n_rows;
    int p=X.n_cols;
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

vec Gate::getZeta(mat z, mat pi){
    int r=pi.n_cols;
    int n=pi.n_rows;
    vec zeta(n*r);
    for(int i=0; i<n; i++){
        mat C=this->makeAchol(vectorise(pi.row(i)));
        mat helper=z.row(i)-pi.row(i);
        helper.reshape(helper.n_cols,helper.n_rows);
        zeta(span(i*r,(i+1)*r-1))=solve(C.t(),helper);
}
return zeta;
}

vec Gate::findGammaQR(mat X, mat z, mat Omega, mat* R){
    int r=z.n_cols;
    int n=X.n_rows;
    int p=X.n_cols;
    mat Omega_chol=chol(Omega);
    vec gamma(p*r);
    gamma.zeros();
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

vec Gate::findGammaQR(mat X, mat z, mat Omega){
    mat R;
    return this->findGammaQR(X,z,Omega,&R);
}

vec Gate::proposeGamma(vec gammaold, mat X, mat z, mat Omega){
   vec mu_gamma(gammaold.size());
   mu_gamma.zeros();
   mat R;
   vec gammahat = this->findGammaQR(X, z, Omega,&R);
   vec v(gammaold.size(),fill::randn);
   mat Sigma=(R.t()*R).i();
   mat Sigma_gamma=Omega.i();
   mat RHS(R.n_rows,1);
   RHS.zeros();
   RHS.rows(0,gammahat.size()-1)=v;
   vec gammanew=gammahat+solve(sqrt(SigmaMultiple)*R,RHS);
   double loglik_old=this->loglik(z,this->pi_calculator(X,gammaold));
   double loglik_new=this->loglik(z,this->pi_calculator(X,gammanew));
   double proposal_old=sum(this->logmvndensity(gammaold,gammahat,Sigma));
   double proposal_new=sum(this->logmvndensity(gammanew,gammahat,Sigma));
   double prior_old=sum(this->logmvndensity(gammaold,mu_gamma, Sigma_gamma));
   double prior_new=sum(this->logmvndensity(gammanew,mu_gamma, Sigma_gamma));
   double acceptance=loglik_new-loglik_old+proposal_old-proposal_new+prior_new-prior_old;
   double u=randu();
   bool accept=u<exp(acceptance);
   cout<<"Accept/Reject:"<<accept<<endl;
  if(accept==1) return gammanew;
  if(accept==0) return gammaold; 
}


vec Gate::logmvndensity(vec response, vec mean, mat Sigma){
   int k = Sigma.n_cols;
   //return 1/(pow(2*M_PI,k/2)*sqrt(det(Sigma)))*exp(-0.5*(response-mean).t()*Sigma.i()*(response-mean)); - not log scale
   return -k/2*log(2*M_PI)-0.5*log(det(Sigma))-0.5*(response-mean).t()*Sigma.i()*(response-mean);
}

int Gate::getDescendantIndex(Node* node){
    //cout<<node->name<<endl;
    for (int i=0;i<this->countChildren();i++){
        if(node==this->Children[i])
            return i;
    }
    Gate* Parent=node->getParent();
    if(Parent==NULL) return -1;
    cout<<Parent->name<<endl;
    //return -99;
   return this->getDescendantIndex(Parent);
}

int Gate::issueIDLR(int start){
    this->idLR=start++;
    for(int i=0;i<countChildren();i++){
        start=this->Children[i]->issueIDLR(start);
    }
    return start;
}

int Gate::issueIDLR(){
    int start=0;
    return this->issueIDLR(start);
}

  
 int Gate::rightMostNodeID(){
     return this->Children[countChildren()-1]->rightMostNodeID();
 }


// vec Gate::getDescendantRange(){
//     vec result;
//     result<<this->leftMostNodeID()<<this->rightMostNodeID();
//     return result;
// }

int Gate::isInRange(Node* node){
    vec range=this->getDescendantRange(); 
    return range[0] <= node->idLR && node->idLR <= range[1];
}

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
