#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <cmath>

using namespace std;

//Global constants
const int D {5};
const int spin_deg {2};
const int mps_size {3};

//necessary structures
struct index_t{
  int dims[3];
};

struct Tensor_t{
  double data[D][D][spin_deg];
  struct index_t index;
};

struct MPS_t{
  struct Tensor_t *sites;
};

//function declarations
Tensor_t makeTensor(index_t theIndex);
Tensor_t matr_mult(Tensor_t t1, Tensor_t t2, int *config);
double evaluate(int *config, MPS_t mps);
MPS_t adjoint(MPS_t mps);
double overlap(MPS_t mps1, MPS_t mps2);
Tensor_t triple_matr_product(Tensor_t t1, Tensor_t t2, Tensor_t t3);

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
}


Tensor_t makeTensor( index_t theIndex){
//  Tensor(Array(Float64,tuple(theIndex.index_dims[:]...)),theIndex)
//  Tensor(rand(tuple(theIndex.index_dims[:]...)),theIndex)
  struct Tensor_t t;
  t.index =  theIndex;

  //fill tensor with random values between 0-1
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




/*multiply bond matrices of two tensors given a physical configuration
  stores resulting matrix in Tensor t at physical index 0
*/
Tensor_t matr_mult(Tensor_t t1, Tensor_t t2, int *config){

  struct Tensor_t t;
  int sum;
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


