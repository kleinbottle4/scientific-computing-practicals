// -*- origami-fold-style: triple-braces; -*-
#ifndef arrow_matrix_dot_H
#define arrow_matrix_dot_H

#include "arrow-addresser.H"
#include <vector>
#include "realvec.H"

class arrow_matrix {
private:
  arrow_addresser ad;

  // Purely to ease reading
  int owner(int m) const {return ad.get_owner()[m];}
  int neighbour(int m) const {return ad.get_owner()[m];}

  int N;
  std::vector<double> upper;
  std::vector<double> lower;
  std::vector<double> diag;
  int number_of_asymmetries;

public:
  arrow_matrix(unstructured_mesh m)// {{{
    : ad(m)
  {
    N = ad.get_n();
    diag.resize(N, 0.0);
    upper.resize(ad.get_num_off_diag_elts() / 2, 0.0);
    lower.resize(ad.get_num_off_diag_elts() / 2, 0.0);
    number_of_asymmetries = 0;
  }
  // }}}

  arrow_matrix(const arrow_matrix &A) // {{{
  {
    ad = A.ad;
    N = A.N;
    upper = A.upper;
    lower = A.lower;
    diag = A.diag;
    number_of_asymmetries = A.number_of_asymmetries;
  }
  // }}}

  bool symmetric () const// {{{
  {
    return (number_of_asymmetries == 0);
  }
  // }}}

  bool diagonally_dominant () const // {{{
  {
    std::vector<double> off_diagonal_sums(N); // One entry for each row

    // <Find the off-diagonal sum of each row. // {{{
    for (int m = 0; m < ad.get_num_off_diag_elts() / 2; m++) {
      int i = owner(m);
      int j = neighbour(m);
      int A_ij = upper[m];
      int A_ji = lower[m];
      off_diagonal_sums[i] += abs(A_ij);
      off_diagonal_sums[j] += abs(A_ji);
    }
    // }}}

    // <Return false unless we have weak diagonal dominance. // {{{
    for (int i=0; i<N; i++){
      if (diag[i] < off_diagonal_sums[i]){
	return false;
      }
    }
    // }}}

    // <Return true if at least one row has strict diagonal dominance. // {{{
    for (int i=0; i<N; i++){
      assert(diag[i] >= off_diagonal_sums[i]);
      if (diag[i] > off_diagonal_sums[i]){
	return true;
      }
    }
    // }}}

    return false;
  }
  // }}}

  // spd{{{
  // Helper functions {{{

  bool diagonals_are_positive() const
  {
    // Check that all diagonal coefficients are positive
    for (int i=0; i<N; i++)
      if (diag[i] <= 0)
        return false;
    return true;
  }

  bool spd_inplace_via_ib_linear_algebra_method()
  {
    for (int i = 0; i < N; i++) {
      // Eliminate all the off-diagonal coefficients in row $i$ {{{
      for (int j = i + 1; j < N; j++) {
	// If [[diag[i] ]] is zero then $e_i^T A e_i = 0$ so $A$ is
	// not positive definite. {{{
	if (diag[i] == 0) {
	  return false;
	}
	// }}}

        // Perform a column operation to eliminate $A_{i, j}$ {{{
        double c = (*this)(i, j) / diag[i];
	// Add c times column i to column j.
        for (int t = 0; t < N; t++)
	  (*this)(t, j) += c * (*this)(t, i);
        // }}}

        // Perform the equivalent row operation. {{{

	// I.e. add c times row i to row j. By the associativity of
	// matrices, the order of the row and column additions is
	// irrelevant.

        for (int t = 0; t < N; t++)
	  (*this)(j, t) += c * (*this)(i, t);
	
        // }}}
      }
      // }}}
    }

    return diagonals_are_positive();

  }
  // }}}

  bool spd() const
  {
    if (!symmetric()){
      return false;
    } else if (diagonally_dominant() && diagonals_are_positive()) {
      return true;
      // By the Gershgorin theorem about the boundaries of the discs.
    } else {
      arrow_matrix copy(*this);
      return copy.spd_inplace_via_ib_linear_algebra_method();
    }
  }
  // }}}


  // Multiply a vector {{{

  struct dimension_mismatch_error {};

  realvec operator*(realvec x)
  {
    // Check dimensions.
    if (x.n != N) {
      throw dimension_mismatch_error();
    }
    // Initialize return vector, b.
    realvec b(N, 0.0);

    // Multiply by iterating through the coefficients of upper and
    // lower simultaneously.

    // Think of the formula b_i = sum{j=1}^n A_{ij} x_j for |upper|

    // and b_j = sum{i=1}^n A_{ji} x_i for |lower|,

    // but actually we are doing b_s = A{ss}*x_s + sum_{t<n/2} A_{ts}
    // + sum_{t>n/2} A_{ts} with (s, t) being (i, j) for upper and (s,
    // t) being (j, i) for lower

    for (int s = 0; s < N; s++)
      b[s] = diag[s];

    for (int m = 0; m < upper.size(); m++){
      int i = owner(m);
      int j = neighbour(m);
      double A_ij = upper[m];
      double A_ji = lower[m];
      b[i] += A_ij * x[j];
      b[j] += A_ji * x[i];
    }

    return b;
  }
  // }}}

  // Resetting and setting coefficients {{{

  // Helper functions for operator() {{{

  // Functions to reset values of existing coefficients {{{
  void reset_upper_coef(int m, double new_value)
  {
    double old_value = upper[m];
    upper[m] = new_value;
    // update [[number_of_asymmetries]]{{{
    if (old_value == lower[m] && new_value != old_value) {
      number_of_asymmetries--;
    } else if (old_value != lower[m] && new_value == old_value) {
      number_of_asymmetries++;
    }// }}}

  }

  void reset_lower_coef(int m, double new_value)
  {
    double old_value = lower[m];
    lower[m] = new_value;
    // update [[number_of_asymmetries]]{{{
    if (old_value == upper[m] && new_value != old_value) {
      number_of_asymmetries--;
    } else if (old_value != lower[m] && new_value == old_value) {
      number_of_asymmetries++;
    }// }}}

  }

  void reset_diag_coef(int m, double new_value)
  {
    diag.at(m) = new_value;
  }
  // }}}

  // Setting a new value {{{
  void set_new_upper_coef(int i, int j, double value)
  {
    ad.fill_in(i, j);
    upper.push_back(value);
    lower.push_back(0.0);
    // update number of assymmetries {{{
    if (value != 0.0){
      number_of_asymmetries++;
    }
    // }}}
  }

  void set_new_lower_coef(int i, int j, double value)
  {
    ad.fill_in(i, j);
    lower.push_back(value);
    upper.push_back(0.0);
    // update number of assymmetries {{{
    if (value != 0.0){
      number_of_asymmetries++;
    }
    // }}}
  }
  // }}}
  // }}}

  double operator()(int i, int j) const // {{{
  {
    // Get A_{i,j}. If not stored then return 0.0. 
    if (i == j) { 
      return diag[i];      
    } else if (ad.stored(i, j)) {
      int m = ad.index(i, j);
      if (j < i){
	return lower[m];
      } else {
	return upper[m];
      }
    } else {
      return 0.0;
    }
  }
  // }}}

  double & operator() (int i, int j) // {{{
  {
    if (i == j){
      return diag[i]; // reference
    } else if (ad.stored(i, j)) {
      int m = ad.index(i, j);
      if (j < i)
	return lower[m]; // reference
      else
	return upper[m]; // reference
    } else if (j < i) {
      set_new_lower_coef(i, j, 0.0);
      return lower[lower.size()-1]; // reference
    } else {
      assert(i < j);
      set_new_upper_coef(i, j, 0.0);
      return upper[upper.size() -1]; // reference
    }
  }
  // }}}

  // Do we need to be careful if we write something like
  // A(1, 2) = A(3, 4) = 5.6;
  // ?
  // }}}

  // Getters for the private data // {{{
  const arrow_addresser & get_ad() const { return ad; }
  const std::vector<double> & get_upper() const { return upper; }
  const std::vector<double> & get_lower() const { return lower; }
  const std::vector<double> & get_diag() const { return diag; }
  int get_number_of_asymmetries() const { return number_of_asymmetries; }
  int dim() const { return N; }
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

#endif
