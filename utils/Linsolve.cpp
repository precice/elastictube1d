
#include "Linsolve.hpp"
#include "AD.hpp"
#include <cmath>
#include <assert.h>

using namespace AD;

namespace LA{
  
void linsolve(
    const int n,
    const int m,
    AD::dualReal** mat,
    AD::dualReal*  rhs,
    AD::dualReal*  x
){
    
  dualReal** Q, **R, *b;
  
  b   = (dualReal*)calloc((m), sizeof(dualReal));    if(b == NULL) exit(0);
  Q = (dualReal**)calloc((n), sizeof(dualReal*));    if(Q == NULL) exit(0);
  for (int i = 0; i < n; ++i){  
    Q[i] = (dualReal*)calloc((m), sizeof(dualReal)); if(Q[i] == NULL) exit(0);
  }
  R = (dualReal**)calloc((m), sizeof(dualReal*));    if(R == NULL) exit(0);
  for (int i = 0; i < m; ++i){
    R[i] = (dualReal*)calloc((m), sizeof(dualReal)); if(R[i] == NULL) exit(0);
  }
  
  // decompose A = QR
  modifiedGramSchmidt(n, m, mat, Q, R);
  
  // compute b = Q^T*rhs
  dualReal tmp;
  for(int j=0; j<m; j++){
    tmp = dualReal(0.);
     for(int i=0; i<n; i++){
      tmp = tmp + Q[i][j] * rhs[i];
    }
    b[j] = tmp;
  }
  
  // backsubst. Rx = b
  backSubstitution(n, m, R, b, x);
  
  // free memory
  free(b);
  for (int i = 0; i < n; ++i)
    free(Q[i]);
  free(Q);
  for (int i = 0; i < m; ++i)
    free(R[i]);
  free(R);
}

void modifiedGramSchmidt (
  const int n,
  const int m,
  dualReal** A,
  dualReal** Q,
  dualReal** R
) {  
  for (int i=0; i < n; i++){
    Q[i][0] = A[i][0];
  }

  for (int k=0; k < m; k++){
    // Compute norm of k-th column vector of (modified) A
    R[k][k] = dualReal(0.);
    for (int i=0; i < n; i++){
      R[k][k] = R[k][k] + AD::pow(A[i][k], 2);
    }
    R[k][k] = AD::sqrt(R[k][k]);

    // Normalize k-th column of matrix Q
    for (int i=0; i < n; i++){
      Q[i][k] = A[i][k] / R[k][k];
    }

    // Compute entries in R and next orthonormal vector
    for (int j=k+1; j < m; j++){
      // Compute entries of R from current Q and A
      for (int i=0; i < n; i++){
	R[k][j] = R[k][j] + Q[i][k] * A[i][j];
      }
      // Subtract contributions from computed to open orthonormal vectors
      for (int i=0; i < n; i++){
	A[i][j] = A[i][j] - Q[i][k] * R[k][j];
      }
    }
  }
}


void backSubstitution (
  const int n,
  const int m,
  dualReal** matrix,
  dualReal*  rhs,
  dualReal*        x
) {
  for (int i = n-1; i >= 0; i--){
    x[i] = rhs[i];
    for (int j = i+1; j < n; j++ ){
      x[i] = x[i] - matrix[i][j] * x[j];
    }
    x[i] = x[i] / matrix[i][i];
  }
}

}