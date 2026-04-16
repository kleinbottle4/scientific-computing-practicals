// -*- origami-fold-style: triple-braces; -*-
#include <vector>
#include "vec.h"

typedef compressed_row_matrix sparse_matrix;

class Axeqb {
  sparse_matrix A;
  vec x;
  vec b;
  vec r;
  int N;
  
  vec getxref()
  {
    double avg = float(x.sum()) / N;
    return vec(N, avg);
  }

public:
  Axeqb(sparse_matrix A, vec x, vec b)
  {
    this->A = A;
    this->x = x;
    this->b = b;
    this->r = b - A*x;
    this->N = A.shape[0];
    // Assert that A is a square matrix
    assert(A.shape[1] == self.N);
    // Assert that x, r and b have length N {{{
    assert(x.size() == N);
    assert(b.size() == N);
    assert(r.size() == N);
    // }}}
  }

  vec getresidual()
  {
    return r;
  }

  void visualizeresidual(double scale)
  {
    for (int j = 0; j < N; j++) {
      for (int a = 0; a < ; a++) {
	std::cout << '*';
      }
      std::cout << '\n';
    }
  }
    
  double getn() const
  {
    // Scalar normalisation factor
    double n = 0.0;
    vec v1 = A*x - A*getxref();
    vec v2 = b - A*getxref();
    for (int i=0; i<N; i++) {
      n += fabs(v1[i]) - fabs(v2[i]);
    }
    return n;
    // What is the physical meaning?
  }

  double onenorm()
  {
    double a = 0.0;
    for (int j = 0; j < N; j++) {
      a += fabs(r[j]);
    }
    return a;
  }

  double twonorm()
  {
    double a = 0.0;
    for (int j = 0; j < N; j++) {
      a += fabs(r[j]) * fabs(r[j]);
    }
    return sqrt(a);
  }

  double infinitynorm()
  {
    double a = 0.0;
    for (int j = 0; j < N; j++) {
      if (r[j] > a) {
	a = r[j];
      }
    }
    return a;
  }

  double getnormalisedscalarresidual() const
  {
    return r / getn();
  }
};
