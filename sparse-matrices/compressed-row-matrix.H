// -*- origami-fold-style: triple-braces; -*-
#ifndef compressed_row_matrix_dot_H
#define compressed_row_matrix_dot_H

#include "realvec.H"
#include "unstructured-mesh.H"
#include "compressed-row-addresser.H"
#include "iostream"

class compressed_row_matrix {
  compressed_row_addresser ad;
  std::vector<double> coeffs;
  int N;
  
  double col(int m) const { return ad.get_col()[m]; }
  double row(int m) const { return ad.get_row()[m]; }

public:

  compressed_row_matrix(unstructured_mesh m) // {{{
    : ad(m)
  {
    N = ad.N();
    coeffs.resize(ad.num_off_diag_elts() + ad.N(), 0.0);
  }
  // }}}

  // Accessors for element (i, j){{{
  double operator()(int i, int j) const
  {
    std::pair<bool, int> a = ad.get_stored_and_index(i, j);
    bool stored = a.first;
    int m = a.second;
    if (stored)
      return coeffs[m];
    else 
      return 0.0;
  }

  struct unauthorized_fill_in {
    int i;
    int j;
  };

  // With fill-in
  double & operator()(int i, int j, bool do_fill_in = false)
  {
    std::pair<bool, int> a = ad.get_stored_and_index(i, j);
    bool stored = a.first;
    int m = a.second;
    if (stored) {
      return coeffs[m];
    } else {
      // Fill-in

      if (!do_fill_in)
	throw unauthorized_fill_in();
      
      // Fill-in for (i, j) in the addresser
      int m = ad.fill_in(i, j);
      // Insert the new (i, j) entry into the correct place (m) in |coeffs|
      coeffs.insert(coeffs.begin() + m, 0.0);
      // Return a reference to the new entry.
      return coeffs[m];
    }
  }
  // }}}

  bool symmetric() const // {{{
  {
    // For the arrow format this would be faster.

    for (int i = 0; i < N; i++) {
      for (int m = row(i); m < row(i + 1); m++) {
        // The opposite coefficient.
        int j = col(m);
        assert(coeffs[m] == (*this)(i, j));
        double b = (*this)(j, i);

        if (coeffs[m] != b) 
          return false;
        
      }
    }
    return true;
  }
  // }}}

  bool diagonally_dominant() const// {{{
  {
    bool at_least_one_strict_inequality = false;

    for (int i = 0; i < N; i++){
      // Return false if row i is diagonally enslaved. {{{
      double mag_diagonal = fabs((*this)(i, i));

      double off_diag_mag_sum = 0.0;

      for (int m = row(i); m < row(i+1); m++)
        off_diag_mag_sum += (col(m) == i ? 0.0 : fabs(coeffs[m]));

      if (mag_diagonal < off_diag_mag_sum)
        return false;
      // }}}

      //  Flag if this row is strictly diagonally dominant. {{{
      if (mag_diagonal > off_diag_mag_sum)
        at_least_one_strict_inequality = true;
      // }}}
    }

    // By this point we know that all the rows are (weakly) diagonally
    // dominant.

    return at_least_one_strict_inequality;
  }
  // }}}

  // spd {{{

  bool diagonals_are_positive() const // {{{
  {
    for (int i = 0; i < N; i++)
      if ((*this)(i, i) <= 0.0)
        return false;

    return true;
  }
  // }}}

  bool spd() const
  {
    if (!symmetric()){
      return false;
    } else if (diagonally_dominant() && diagonals_are_positive()) {
      // By the Gershgorin theorem about the boundaries of the discs.
      return true;
    } else {
      // Long method from IB Linear Algebra {{{

      compressed_row_matrix A = *this;

      for (int i = 0; i < N; i++){

        // Edge case {{{
        if (A(i, i) == 0.0)
          return false;
        // }}}

        // Now we operate on A via A -> B^T A B

        // Clear the i-th row and i-th column, except A_{ii} {{{

        for (int j = i + 1; j < N; j++) {
          // Eliminate A_{ij} {{{
          // I.e. Subtract c * column i from column j.
          double c = A(i, j) / A(i, i);

          for (int t = 0; t < N; t++)
            A(t, j) -= c * A(t, i);

          // }}}

          // Perform the transpose operation {{{
          //   I.e. Subtract c * row i from row j.
          //   This eliminates A_{ji}.
          for (int t = 0; t < N; t++)
            A(j, t) -= c  *A(t, i);
          // }}}
        }
        // }}}

      }

      // By this point, A has become diagonal.
      // If the diagonal are all positive then A is positive definite.
      return diagonals_are_positive();

      // }}}
    }
  }
  // }}}

  // Multiply with a vector{{{

  struct matrix_and_vector_dimensions_incompatible{};

  realvec operator*(realvec x)
  {
    // Check dimensions{{{
    if (x.n != N) 
      throw matrix_and_vector_dimensions_incompatible();
    // }}}

    realvec b(N, 0.0); // Initially all zeroes.

    for (int i = 0; i < N; i++){
      // Use the usual formula b_i = sum_{j=1}^n A_{ij} x_j
      //  but skip zeroes.
      // {{{
      for (int m = row(i); m < row(i + 1); m++){
        int j = col(m);
        double A_ij = coeffs[m];
        b[i] += A_ij * x[j];
      }
      // }}}
    }

    return b;
  }
  // }}}

  // Getters for private data{{{
  int dim() const { return N; }
  compressed_row_addresser get_ad() const { return ad; }
  std::vector<double> get_coeffs() const { return coeffs; }
  // }}}


  // Printing {{{
  void print() const
  {
    std::cout << "[\n";
    for (int i=0; i < dim(); i++){
      std::cout << "[";
      for (int j=0; j < dim(); j++){
	std::cout << (*this)(i, j) << ", ";
      }
      std::cout << "];\n";
    }
    std::cout << "]\n";
  }
  // }}}

};


#endif // compressed_row_matrix_dot_H
