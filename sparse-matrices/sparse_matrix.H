#include "realvec.H"

class sparse_matrix {
  virtual double operator()(int i, int j) const = 0;
  virtual double & operator()(int i, int j) = 0;
  virtual realvec operator*(realvec) const = 0;
  virtual bool symmetric() const = 0;
  virtual bool spd() const = 0;
  virtual bool dim() const = 0;
};

// Currently I do not use this class, but I obey it in style.
