#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

//Global variables
const int D = {5};
const int Dw = {5};
const int spin_deg = {2};
const int mps_size = {3};

//necessary structures
struct index_t{
  int dims[4];
};

struct Tensor_t{
  double data[D][D][spin_deg];
  struct index_t index;
};

//tensor of product before compression
struct pTensor_t{
  double data[D*Dw][D*Dw][spin_deg];
  struct index_t index;
};

//Tensor in an MPO
struct oTensor_t{
  double data[Dw][Dw][spin_deg][spin_deg];
  struct index_t index;
};

struct MPS_t{
  struct Tensor_t *sites;
};

//mps resulting from product before compression
struct pMPS_t{
  struct pTensor_t *sites;
};

struct MPO_t{
  struct oTensor_t *sites;
};

//function declarations
Tensor_t makeTensor(index_t theIndex);
MPO_t makeMPO();
Tensor_t matr_mult(Tensor_t t1, Tensor_t t2, int *config);
double evaluate(int *config, MPS_t mps);
MPS_t adjoint(MPS_t mps);
double overlap(MPS_t mps1, MPS_t mps2);
Tensor_t triple_matr_product(Tensor_t t1, Tensor_t t2, Tensor_t t3);
MPS_t apply_one_site_H(int H[spin_deg][spin_deg], int i, MPS_t mps);
pMPS_t apply_MPO(MPO_t o, MPS_t m);
int base_n_counter(int *count);

int main(){

  struct MPS_t mps1, mps2;
  mps1.sites = new Tensor_t[mps_size];
  mps2.sites = new Tensor_t[mps_size];

  //initialize first Tensor
  struct index_t ind;
  ind.dims[0] = 1;
  ind.dims[1] = D;
  ind.dims[2] = spin_deg;
  mps1.sites[0] = makeTensor(ind);
  mps2.sites[0] = makeTensor(ind);

  //initialize last Tensor
  ind.dims[0] = D;
  ind.dims[1] = 1;
  mps1.sites[mps_size-1] = makeTensor(ind);
  mps2.sites[mps_size-1] = makeTensor(ind);

  ind.dims[1] = D;

  //initialize all other tensors in the mps
  for(int i=1; i<mps_size-1;i++){
    mps1.sites[i] = makeTensor(ind);
    mps2.sites[i] = makeTensor(ind);
  }

  printf("MPS generated\n");
  
  int config[3] = {0,1,1};
  double eval1 = evaluate(config, mps1);
  printf("MPS 1 Evaluated: %f\n", eval1);
  
  double eval2 = evaluate(config, mps2);
  printf("MPS 2 Evaluated: %f\n", eval2);

  double ovlp = overlap(mps1, mps2);
  printf("Overlap: %f\n", ovlp);

  MPO_t o = makeMPO();
  printf("MPO Generated\n");

  pMPS_t m = apply_MPO(o, mps1);

  if( m.sites[0].index.dims[0] == -1){
    printf("MPO-MPS Application Failed\n");
  }
  else{
    printf("MPO-MPS Application Successful\n");
  }
}


Tensor_t makeTensor( index_t theIndex){
//  Tensor(Array(Float64,tuple(theIndex.index_dims[:]...)),theIndex)
//  Tensor(rand(tuple(theIndex.index_dims[:]...)),theIndex)
  struct Tensor_t t;
  t.index =  theIndex;

  //fill tensors
  for(int i=0; i<theIndex.dims[0]; i++){
    for(int j=0; j<theIndex.dims[1]; j++){
      for(int k=0; k<theIndex.dims[2];k++){
        //t.data[i][j][k] = double(rand()/RAND_MAX);
        t.data[i][j][k] = 1;
      }
    }
  }
  return t;
}

MPO_t makeMPO(){

  struct MPO_t m;
  m.sites = new oTensor_t[mps_size];

  //loop through sites of the MPO
  for(int i=0; i<mps_size; i++){

    //assign tensor dimensions
    m.sites[i].index.dims[0] = Dw;
    m.sites[i].index.dims[1] = Dw;
    m.sites[i].index.dims[2] = spin_deg;
    m.sites[i].index.dims[3] = spin_deg;

    //fill tensor
    for(int j=0; j<m.sites[i].index.dims[0]; j++){
      for(int k=0; k<m.sites[i].index.dims[1]; k++){
        for(int l=0; l<m.sites[i].index.dims[2];l++){
          for(int n=0; n<m.sites[i].index.dims[3];n++){
            m.sites[i].data[j][k][l][n] = 1;
          }
        }
      }
    }
  }
  return m;
}



/*multiply bond matrices of two tensors given a physical configuration
  stores resulting matrix in Tensor t at physical index 0
*/
Tensor_t matr_mult(Tensor_t t1, Tensor_t t2, int *config){

  int sum;

  struct Tensor_t t;
  t.index.dims[0] = t1.index.dims[0];
  t.index.dims[1] = t2.index.dims[1];
  t.index.dims[2] = 1;

  //loop through every element in resulting matrix
  for(int i = 0; i < t.index.dims[0]; i++){
    for(int j = 0; j < t.index.dims[1]; j++){

      //calculate each element of the matrix with matrix multiplication
      sum = 0;
      for(int k = 0; k < t1.index.dims[1]; k++){
        sum += t1.data[i][k][config[0]]*t2.data[k][j][config[1]];
      }
      t.data[i][j][0] = sum;
    }
  }
  return t;
}


/*Evaluates an mps given a physical configuration
  the configuration is 1-indexed
*/
double evaluate(int *config, MPS_t mps){

  struct Tensor_t t;

  for(int i=0; i < mps_size-1; i++){

    if(i == 0){
      t = matr_mult(mps.sites[0], mps.sites[1], config);
    }
    else{
      t = matr_mult(t, mps.sites[i+1], &config[i]);
    }
    config[i+1] = 0;
  }

  //sanity check
  if(t.index.dims[0] == t.index.dims[1] && t.index.dims[0]== 1){

    //final product will be stored here
    return t.data[0][0][0];
  }
  return -1;
}

//takes adjoint of every bond matrix at every site of an mps
//currently assummes every value is real i.e. it takes the transpose
MPS_t adjoint(MPS_t mps){

  struct MPS_t ret;
  ret.sites = new Tensor_t[mps_size];

  //loop through each site of mps
  for(int i=0; i<mps_size; i++){

    ret.sites[i].index.dims[0] = mps.sites[i].index.dims[1];
    ret.sites[i].index.dims[1] = mps.sites[i].index.dims[0];
    ret.sites[i].index.dims[2] = spin_deg;

    //loop through all physical indices
    for(int j = 0; j < spin_deg; j++){

      //case if matrix is not a column vector
      if(ret.sites[i].index.dims[0] <= ret.sites[i].index.dims[1]){

        //loop through upper diagonal elements of the bond matrix
        for(int k=0; k<ret.sites[i].index.dims[0]; k++){
          for(int l = k; l<ret.sites[i].index.dims[1]; l++){

            //copy elements from accross the diagona
            ret.sites[i].data[k][l][j] = mps.sites[i].data[l][k][j];
            ret.sites[i].data[l][k][j] = mps.sites[i].data[k][l][j];
          }
        }
      }
      //case if matrix is a column vector
      else{

        //loop through first column
        for(int k = 0; k<ret.sites[i].index.dims[1]; k++){
          ret.sites[i].data[0][k][j] = mps.sites[i].data[k][0][j];
        }
      }
    }
  }
  return ret;
}

double overlap(MPS_t mps1, MPS_t mps2){

  //stores current values in computation
  struct Tensor_t t;

  struct MPS_t adj = adjoint(mps2);

  //create identity matrix
  //only needs to be 1 dimensional
  struct Tensor_t id;
  id.index.dims[0] = 1;
  id.index.dims[1] = 1;
  id.index.dims[2] = 1;
  id.data[0][0][0] = 1;
  

  //loop over all sites
  for(int i=0; i<mps_size; i++){

    //have identity as middle term if first instance of loop
    if(i == 0){
      t = triple_matr_product(adj.sites[i], id, mps1.sites[i]);
    }
    //otherwise feed in result from previous loop
    else{
      t = triple_matr_product(adj.sites[i], t, mps1.sites[i]);
    }

    //check computation was successful
    if(t.index.dims[0] == -1){
      return -1;
    }

    //sum over all physical indices, store result in t at physical index 0
    if(spin_deg > 1){
      for(int j=1; j<spin_deg; j++){
        for(int k=0; k<t.index.dims[0]; k++){
          for(int l=0; l<t.index.dims[1]; l++){
            t.data[k][l][0] += t.data[k][l][j];
          }
        }
      }
    }
  }

  if(t.index.dims[0] != 1 || t.index.dims[1] != 1){
    return -1;
  }

  return t.data[0][0][0];
}

/*
  Computes the product of three bond matrices over all physical indices
  Stores result in physical index's corresponding bond matrix in t
  Middle tensor only evaluated at physical index 0, product of previous calculation or identity
*/
Tensor_t triple_matr_product(Tensor_t t1, Tensor_t t2, Tensor_t t3){

  Tensor_t t, temp;
  t.index.dims[0] = t1.index.dims[0];
  t.index.dims[1] = t3.index.dims[1];
  t.index.dims[2] = spin_deg;

  //ensure correct dimensions for matrix multiplication
  if(t1.index.dims[1] != t2.index.dims[0] || t2.index.dims[1] != t3.index.dims[0]){
    t.index.dims[0] = -1;
    return t;
  }

  int config[2] = {0,0};

  //loop over all physical indices
  for(int i=0; i<spin_deg; i++){

    //only evaluate t1, t3 at physical index
    config[0] = i;
    config[1] = 0;
    temp = matr_mult(t1, t2, config);
    config[0] = 0;
    config[1] = i;
    temp = matr_mult(temp, t3, config);

    //store result in corresponding bond matrix in t
    for(int j=0; j<t.index.dims[0]; j++){
      for(int k=0; k<t.index.dims[1]; k++){
        t.data[j][k][i] = temp.data[j][k][0];
      }
    }
  }
  return t;
}
/*
  Apply H(spin_deg X spin_deg matrix) to site i(0 indexed) of the MPS
*/
MPS_t apply_one_site_H(int H[spin_deg][spin_deg], int i, MPS_t mps){

  MPS_t ret;
  ret.sites = new Tensor_t[mps_size];

  //loop through all sites
  for(int j=0; j<mps_size; j++){

    //new mps has same dimensions as old one
    ret.sites[j].index.dims[0] = mps.sites[j].index.dims[0];
    ret.sites[j].index.dims[1] = mps.sites[j].index.dims[1];
    ret.sites[j].index.dims[2] = mps.sites[j].index.dims[2];

    //loop through bond matrix indices
    for(int k=0; k<ret.sites[j].index.dims[0]; k++){
      for(int l=0; l<ret.sites[j].index.dims[1]; l++){

        //if not applying H to this site, just copy values
        if(j != i){
          for(int n=0; n<ret.sites[j].index.dims[2]; n++){
            ret.sites[j].data[k][l][n] = mps.sites[j].data[k][l][n];
          }
        }

        //otherwise apply the operator
        else{
          for(int m=0; m<spin_deg; m++){
            for(int n=0; n<spin_deg; n++){
              ret.sites[j].data[k][l][m] += H[m][n]*mps.sites[j].data[k][l][n];
            }
          }
        }
      }
    }
  }
  return ret;
}


pMPS_t apply_MPO(MPO_t o, MPS_t m){

  pMPS_t ret;
  ret.sites = new pTensor_t[mps_size];

  //loop through sites of MPS
  for(int i=0; i<mps_size; i++){

    ret.sites[i].index.dims[0] = o.sites[i].index.dims[0]*m.sites[i].index.dims[0];
    ret.sites[i].index.dims[1] = o.sites[i].index.dims[1]*m.sites[i].index.dims[1];
    ret.sites[i].index.dims[2] = spin_deg;

    //indices of bond matrix j and k (base 10) can be represented as coupled indices of original matrix ind1 and ind1 
    int ind1[2] = {0,0};
    int ind2[2] = {0,0};

    //loop through bond matrix of resulting MPS
    for(int j=0; j<ret.sites[i].index.dims[0]; j++){

      //reset k counter
      ind2[0] = 0;
      ind2[1] = 0;

      for(int k=0; k<ret.sites[i].index.dims[1]; k++){

        //loop through physical index of resulting MPS (sigma)
        for(int l=0; l<spin_deg; l++){

          //zero out value beforehard
          ret.sites[i].data[j][k][l] = 0;

          //loop over physical index being summed over (sigma')
          for(int n=0; n<spin_deg; n++){
            ret.sites[i].data[j][k][l] += o.sites[i].data[ind1[0]][ind2[0]][l][n]*m.sites[i].data[ind1[1]][ind2[1]][n];
          }
        }

        //increment ind 2 counter, perform check it agrees with k
        if( base_n_counter(ind2) != k){
          ret.sites[0].index.dims[0] = -1;
          return ret;
        }
      }

      //increment ind1 counter, perform check it agrees with j
      if( base_n_counter(ind1) != j){
        ret.sites[0].index.dims[0] = -1;
        return ret;
      }
    }
  }
  return ret;
}

/*
  increment the counter of coupled index, index 0 has goes to Dw, index 1 goes to D
  return current count in base 10, returns -1 when reached end of the count
*/
int base_n_counter(int *count){

  int dec = count[0]*Dw + count[1];

  //return -1 if reached the end of the counter
  if(count[0] == Dw){
    return -1;
  }

  //increment the counter
  //higher order digit stored at index 1
  count[1]++;
  if(count[1] == D){
    count[1] = 0;
    count[0]++;
  }
  return dec;
}
