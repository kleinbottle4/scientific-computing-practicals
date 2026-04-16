std::array<double, n> compute_eigenvalues(double tolerance)
{
  // Use power method and gram-schmidt. Stop when the difference
  // between successive guesses is less than 'tolerance'

  // Maybe we could use repeated squaring.  But this would require
  // us to store another copy of the matrix.  And we may converge
  // very quickly onto an eigevector which we have already found.

  // Obviously not returning the eigenvectors is wasteful, but we do
  // not require them at the present time.


  // PROBLEM!!!!
  // How are we meant to store 100,000,000 eigenvectors?

  std::vector<realvec<n>> known_eigenvectors;
  std::vector<double> known_eigenvalues;
  while (known_eigenvalues.size() < n){
    realvec<n> x;
    for (int i=0; i<n; i++)
      x[i] = 1;
    realvec<n> old_x;
    do {
      old_x = x;
      x = (*this) * x;
      x.gram_schmidt(known_eigenvectors);
      x.renormalise();
    } while ((x - old_x).norm() < tolerance);
    known_eigenvectors.push_back(x);
    known_eigenvalues.push_back(((*this)*x).norm());
  }
  return known_eigenvalues;
}
