// -*- origami-fold-style: triple-braces; -*-
#ifndef vec_dot_H
#define vec_dot_H

#include <iostream>
#include <cmath>
#include <vector>

struct vec {
  int n;
  std::vector<double> coeffs;

  struct dimension_mismatch{};

  // Constructors{{{
  vec() {}
  vec(int nn) : n(nn), coeffs(nn) {}
  vec(int nn, double val) : n(nn), coeffs(nn, val) {}
// }}}

  // Accessors
  double operator[](int i) const { return coeffs[i]; }
  double & operator[](int i) { return coeffs[i]; }

  // Arithmetic {{{
  void operator-=(vec y)
  {
    // Check dimensions agree.{{{
    if (y.n != n)
      throw dimension_mismatch(); 
     // }}}
    for (int i = 0; i < n; i++)
      coeffs[i] -= y[i];
  }

  vec operator-(vec y) const {
    vec x = *this;
    x -= y;
    return x;
  }
  // }}}

  // Dot product {{{
  double dot(vec v) const
  {
    // Check dimensions agree.{{{
    if (v.n != n)
      throw dimension_mismatch();
    // }}}

    double a = 0;
    for (int i = 0; i < n; i++)
      a += coeffs[i] * v.coeffs[i];
    return a;
  }

  double operator|(const vec & v) const
  {
    return this->dot(v);
  }

  void project_onto_complement_of(std::vector<vec> S)
  {
    for (vec v : S)
      *this -= (*this | v);
  }
  // }}}

  // Norms {{{
  double norm() const
  {
    // Standard Euclidean norm
    double sum = 0;
    for (int i = 0; i < n; i++) {
      double a = coeffs[i];
      sum += a * a;
    }
    return sqrt(sum);
  }

  struct cannot_normalize_zero {};

  void renormalise()
  {
    double a = norm();
    if (a == 0)
      throw cannot_normalize_zero();
    for (int i = 0; i < n; i++)
      coeffs[i] /= a;
  }
  // }}}

};

// Printing{{{
std::ostream & operator<<(std::ostream & os, vec v)
{
  std::cout << "{";
  for (int i = 0; i < v.n; i++){
    std::cout << v[i] << ", ";
  }
  std::cout << "}";
  return os;
}
// }}}

#endif // vec_dot_H
