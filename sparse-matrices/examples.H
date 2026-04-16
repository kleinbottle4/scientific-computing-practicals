#include "compressed-row-matrix.H"
#include "compressed-row-addresser.H"
#include "unstructured-mesh.H"
#include "realvec.H"
#include <cassert>
#include "arrow-matrix.H"

unstructured_mesh mesh_example_cartesian_square(int n) // {{{
{
  // Create a Cartesian mesh on a square with n * n cells.
  /* E.g., for n = 4, {{{
   * |----+----+----+----|
   * | 12 | 13 | 14 | 15 |
   * |----+----+----+----|
   * |  8 |  9 | 10 | 11 |
   * |----+----+----+----|
   * |  4 |  5 |  6 |  7 |
   * |----+----+----+----|
   * |  0 |  1 |  2 |  3 |
   * |----+----+----+----|
   * }}} */

  // We only consider topological information.

  unstructured_mesh mesh;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      cell c;
      c.index = n*i + j;

      // Find the neighbours of c:{{{
      for (int d : {-1, 1}){
	for (bool up : {true, false}){
	  int new_i = i + (up ? 0 : d);
	  int new_j = j + (up ? d : 0);

	  if (0 <= new_i && new_i < n && 0 <= new_j && new_j < n) {
	    int m = n * new_i + new_j;
	    c.neighbours.push_back(m);
	  }
	}
      }
// }}}

      mesh.cells.push_back(c);
    }
  }
  return mesh;
}
// }}}

std::ostream & operator<<(std::ostream & os, const std::vector<int> & v)
{
  std::cout << "{ ";
  for (int n : v) {
    std::cout << n << " , ";
  }
  std::cout << " }\n";
  return os;
}

std::ostream & operator<<(std::ostream & os, const std::list<int> & v)
{
  std::cout << "{ ";
  for (int n : v) {
    std::cout << n << " <-> ";
  }
  std::cout << " }\n";
  return os;
}

std::ostream & operator<<(std::ostream & os, const unstructured_mesh & m)
{
  std::cout << "Cell" << "\t" << "Neighbours\n";
  for (cell c : m.cells) {
    std::cout << c.index << '\t' << c.neighbours << '\n';
  }
  return os;
}

static int within_range(int n, int i, int j)
{
  return (0 <= i) && (i < n) && (0 <= j) && (j < n);
}

static int f(int n, int i, int j)
{
  assert(within_range(n, i, j));
  return n*i + j;
}

static int on_the_boundary(int n, int i, int j)
{
  return (i == 0) || (i == n-1) || (j == 0) || (j == n-1);
}

// Make a matrix that can later be used to solve the Poisson equation
// This currently does not deal properly with boundary conditions.
// {{{

compressed_row_matrix matrix_example_poisson(int n)
{
  unstructured_mesh mesh = mesh_example_cartesian_square(n);
  compressed_row_matrix A(mesh);

  for (cell c : mesh.cells) {
    int i = c.index;
    A(i, i) = -4;
    for (int j : c.neighbours) {
      A(i, j) = 1;
    }
  }
  return A;
}

arrow_matrix arrow_matrix_example_poisson(int n)
{
  unstructured_mesh mesh = mesh_example_cartesian_square(n);
  arrow_matrix A(mesh);

  for (cell c : mesh.cells) {
    int i = c.index;
    A(i, i) = -4;
    for (int j : c.neighbours) {
      A(i, j) = 1;
    }
  }
  return A;
}

// }}}
