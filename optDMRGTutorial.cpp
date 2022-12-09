#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

using namespace std;
using Eigen::Tensor;
using Eigen::MatrixXd;
using Eigen::VectorXd;

//Global variables
const int D = {5};
const int Dw = {5};
const int spin_deg = {2};
const int mps_size = {3};

//function declarations
double evaluate(Tensor<double,4> mps, int *config);
MatrixXd get_bond_matr(Tensor<double,4> mps, int site, int phys_index);
Tensor<double,4> t_adjoint(Tensor<double,4> mps);
double overlap(Tensor<double,4> mps2, Tensor<double,4> mps);
Tensor<double,4> apply_MPO(Tensor<double,4> mps, Tensor<double,5> mpo);
Tensor<double,4> compress(Tensor<double,4> mps);

int main(){
  //store mps as degree 4 tensor, first index for the site, next two are bond indices, last is physical index
  Tensor<double,4> mps1(mps_size,D,D,spin_deg);
  mps1.setRandom();

  Tensor<double,4> mps2(mps_size,D,D,spin_deg);
  mps2.setRandom();

  //make first/last matrices a row/column vector
  for(int i=1; i < D; i++){
    for(int j=0; j < D; j++){
      for(int k=0; k < spin_deg; k++){
        mps1(0, i, j, k) = 0;
        mps1(mps_size-1, j, i, k) = 0;
        mps2(0, i, j, k) = 0;
        mps2(mps_size-1, j, i, k) = 0;
      }
    }
  }

  printf("MPS generated successfully\n");

  int config[3] = {0,1,1};
  double eval = evaluate(mps1, config);
  printf("MPS evaluated at [0,1,1]: %f\n",eval);

  double o = overlap(mps1, mps2);
  printf("Overlap performed: %f\n",o);

}

double evaluate(Tensor<double,4> mps, int *config){

  //store product in this matrix, set equal to the identity
  MatrixXd t(D,D);
  t.setZero();
  for(int i=0; i<mps_size; i++){
    t(i,i) = 1;
  }

  //stores bond matrix being multiplied
  MatrixXd b(D,D);

  for(int i=0; i<mps_size; i++){

    b = get_bond_matr(mps, i, config[i]);
    //multiply bond matrix with current product
    t = t*b;
  }

  return t(0,0);
}

/*
  returns mps bond matrix for specific site and physical index
*/
MatrixXd get_bond_matr(Tensor<double,4> mps, int site, int phys_index){

  MatrixXd t(D,D);

  for(int i=0; i < D; i++){
    for(int j = 0; j < D; j++){
      t(i,j) = mps(site, i, j, phys_index);
    }
  }

  return t;
}

/*
  returns mps with adjoint taken on all bond matrices
  first/last sites will switch from row/column vector to column/row vector
*/
Tensor<double,4> t_adjoint(Tensor<double,4> mps){

  MatrixXd b(D,D);
  Tensor<double,4> t(mps_size,D,D,spin_deg);

  //loop through sites and physical index
  for(int i=0; i<mps_size; i++){
    for(int j=0; j<spin_deg; j++){

      //get current bond matrix in loop, take adjoint
      b = get_bond_matr(mps, i, j);
      b.adjointInPlace();

      //copy matrix into return tensor
      for(int m=0; m<D; m++){
        for(int n=0; n<D; n++){
          t(i,m,n,j) = b(m,n);
        }
      }
    }
  }
  return t;
}

/*
  Returns the inner product of two mps
  Adjoint is taken on the first argument
*/
double overlap(Tensor<double,4> mps2, Tensor<double,4> mps){

  Tensor<double,4> adj(mps_size,D,D,spin_deg);
  adj = t_adjoint(mps2);

  //structs to store matrices for multiplication
  MatrixXd tot(D,D);
  MatrixXd b1(D,D);
  MatrixXd b2(D,D);
  MatrixXd pr(D,D);

  //stores results of multiplication before summation
  Tensor<double,3> temp(D,D,spin_deg);

  //loop through all the sites
  for(int i=0; i<mps_size; i++){

    //perform summation over all physical indices
    //this loop fills temp struct with values from multiplication
    for(int j=0; j<spin_deg; j++){

      //create matrix from first sites (column times row vector)
      if(i == 0){

        for(int m=0; m<D; m++){
          for(int n=0; n<D; n++){
            temp(m,n,j) = mps(i,0,m,j)*adj(i,n,0,j);
          }
        }
      }

      //final multiplication (row times column vector) results in scalar
      //store values in (0,0) of bond matrices in temp
      else if(i == mps_size-1){
        b1 = get_bond_matr(adj,i,j);
        b2 = get_bond_matr(mps,i,j);

        //create vectors from the mps structs
        VectorXd r(D);
        VectorXd c(D);
        for(int k=0; k<D; k++){
          r(k) = adj(i,0,k,j);
          c(k) = mps(i,k,0,j);
        }

        //perform multiplication
        VectorXd p(D);
        p = tot*c;
        temp(0,0,j) = r.adjoint()*p;
      }

      //otherwise just perform multiplication with previous result and current bond matrices
      else{
        b1 = get_bond_matr(adj,i,j);   
        b2 = get_bond_matr(mps,i,j);
        pr = b1*tot*b2;

        //copy values into temp struct
        for(int m=0; m<D; m++){
          for(int n=0; n<D; n++){
            temp(m,n,j) = pr(m,n);
          }
        }
      }
    }

    //sum matrices from all physical indeces to update tot
    tot.Zero(D,D);
    for(int j=0; j<spin_deg; j++){
      for(int m=0; m<D; m++){
        for(int n=0; n<D; n++){
          tot(m,n) += temp(m,n,j);
        }
      }
    }
  }
  //return final sum
  return tot(0,0);
}
/*
  returns resulting mps from application of an mpo
  performs the tensor product as well as the renornalization back to
  an mps with bond dimension D
*/
Tensor<double,4> apply_MPO(Tensor<double,4> mps, Tensor<double,5> mpo){

  //stores uncompressed product
  Tensor<double,4> prod(mps_size,D*Dw,D*Dw,spin_deg);

  //struct to store results of tensor contraction at each block
  Tensor<double,3> temp(D,D,spin_deg);

  Tensor<double,2> W(spin_deg,spin_deg);
  Tensor<double,3> M(D,D,spin_deg);

  Tensor<double,4> w4(D,D,spin_deg,spin_deg);
  Tensor<double,3> w3(D,spin_deg,spin_deg);



  //loop through sites of mps
  for(int i=0; i<mps_size; i++){

    //fill blocks of larger matrix
    //these loops are over global indices of block location
    for(int m=0; m<Dw; m++){
      for(int n=0; n<Dw; n++){

        //create tensors to perform contraction
        M = mps.chip(i,0);
        w4 = mpo.chip(i,0);
        w3 = w4.chip(m,0);
        W = w3.chip(n,0);

        Eigen::array<Eigen::IndexPair<int>,1> dims = {Eigen::IndexPair<int>(2,1)};
        temp = M.contract(W, dims);

        //fill block in larger tensor with result of contraction
        for(int l=0; l<spin_deg; l++){
          for(int j=0; j<D; j++){
            for(int k=0; k<D; k++){
              prod(i,m*D+j,n*D+k,l) = temp(j,k,l);     
            }
          }
        }
      }
    }
  }


  Tensor<double,4> ret(mps_size,D,D,spin_deg);
  ret = compress(prod);
  return ret;   
}

/*
  Compresses MPS with bond dimension D*Dw down to dimension D
  Done by performing SVD
*/
Tensor<double,4> compress(Tensor<double,4> mps){

  Tensor<double,4> ret(mps_size,D,D,spin_deg);
  Tensor<double,3> site(D,D,spin_deg);

  for(int i=0; i<mps_size; i++){
    
  }

  return ret;
}