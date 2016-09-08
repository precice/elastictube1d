#ifndef LINSOLVE_H_
#define LINSOLVE_H_

#include "AD.hpp"
#include <boost/concept_check.hpp>

namespace LA{

  /**
   * Produces an economy-size QR decomposition of a matrix A, A is changed.
   *
   * Preconditions:
   * - The column vectors of A have to be linear independent. Otherwise, see
   *   postconditions "R".
   * - If A has n columns and m rows, then:
   *   * n <= m for A
   *   * Q has to have same size than A
   *   * R has to be a n * n matrix
   *
   * Postconditions:
   * - A is changed
   * - Q holds an orthonormal basis of range(A) in its column vectors.
   * - R is an upper triangular matrix, with all main diagonal element unequal
   *   to zero iff A has full rank.
   *
   * The algorithm is taken from the book Numerical Mathematics, Second Edition,
   * by Alfio Querteroni, Riccardo Sacco, and Fausto Saleri, Springer, 2007,
   * pages 85-87.
   */
  
  void linsolve(
    const int n,
    const int m,
    AD::dualReal** mat,
    AD::dualReal*  rhs,
    AD::dualReal*  x
   );

  void modifiedGramSchmidt (
    const int n,
    const int m,
    AD::dualReal** A,
    AD::dualReal** Q,
    AD::dualReal** R
  );
    
  void backSubstitution (
    const int n,
    const int m,
    AD::dualReal** matrix,
    AD::dualReal*  rhs,
    AD::dualReal*        x
  );

}

//#include "GramSchmidt.cpph"

#endif /* LINSOLVE_HPP */
