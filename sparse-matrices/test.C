#include "examples.H"
#include "arrow-matrix.H"

int main()
{
  auto mesh = mesh_example_cartesian_square(4);

  compressed_row_addresser ad(mesh);
  ad.visualize_sparseness_pattern();
  compressed_row_matrix A = matrix_example_poisson(4);
  A.print();

  arrow_addresser bd(mesh);
  bd.visualize_sparseness_pattern();
  arrow_matrix B = arrow_matrix_example_poisson(4);
  B.print();

  // Now let's test an arrow matrix where the sparseness pattern is symmetric
  // but not the coefficients.

  B(1, 2) = -5;
  B.print();
  
  return 0;
}
