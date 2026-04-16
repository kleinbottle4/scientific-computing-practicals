// -*- origami-fold-style: triple-braces; -*-
#ifndef arrow_addresser_dot_H
#define arrow_addresser_dot_H

#include <vector>
#include "unstructured-mesh.H"
#include <cassert>
#include <iostream>

/* Aim {{{
 * Suppose we are given an unstructed mesh and we want to define
 * some matrix with coefficients $A_{i,j}$ for each pair of cells i and
 * j.  This class seeks to store all the connectivity information, for
 * use by the solver class later. This means store the ~owner~ and
 * ~neighbour~ arrays. It also define a function from (i, j) to the index
 * of ~upper~ (or ~lower~ or ~diagonal~) where $A_{i, j}$ will be stored
 * in the matrix class.
 * 
 * I want it to be easy to copy a ~arrow_addresser~.  This is so that we
 * can have many matrices (of the same format) for the same mesh, without
 * having to re-calculate the addressing each time.
 * }}}
 */

class arrow_addresser {
private:
  std::vector<int> owner;
  std::vector<int> neighbour;
  int n; // The matrix is an n by n matrix.
  int num_off_diag_elts;
  int matrix_band;
public:

  // Constructors

  arrow_addresser() {}

  arrow_addresser(const unstructured_mesh &m)// {{{
  {

    /* {{{ The constructor takes an unstructured mesh as input stores the
       connectivity of the mesh in ~owner~ and ~neighbour~. Specifically, if
       mesh cell ~i~ is a neighbour of mesh cell ~j~, then in a matrix we
       would have the $(i, j)$-th and $(j, i)$-th coefficients being
       non-zero. In the addresser, we have ~owner[m]~ set to $i$ and
       ~neighbour[m]~ set to $j$, for some value of $m$, and for another
       value ~mm~, we have ~owner[mm]~ set to $j$ and ~neighbour[mm]~ set to
       $i$. }}} */

    n = m.cells.size();
    num_off_diag_elts = 0;
    matrix_band = 0;
    for (cell c : m.cells) {
      int i = c.index;
      for (int j : c.neighbours) {
	if (i < j) { // upper triangle
	  owner.push_back(i);
	  neighbour.push_back(j);
	}
	if (i != j) 
	  num_off_diag_elts++;
	
	matrix_band = std::max(matrix_band, abs(j - i));
      }
    }
    // Perform a check of the mesh. {{{
    // Since cell i is a neighbour of cell j iff cell j is a neighbour of cell i,
    // we check that the number of off-diagonal elements is even.  If not then
    // there is a problem in the mesh.
    assert(num_off_diag_elts % 2 == 0);  
    // }}}
  }
  // }}}

  // Helper functions for getting the index from (i, j){{{
  bool stored_in_upper(int i, int j) const
  {
    /*
    | %%% |     |     | 3,0 |
    |     | %%% |     |     |
    |     | 2,1 | %%% |     |
    |     |     |     | %%% |
    */
    
    if (j <= i)
      return false;

    for (int m=0; m < num_off_diag_elts / 2; m++)
      if (owner[m] == i && neighbour[m] == j)
	return true;

    return false;
  }

  bool stored_in_lower(int i, int j) const
  {
    // Since the sparseness pattern is symmetric:
    return stored_in_upper(j, i);
  } 

  bool stored(int i, int j) const
  {
    if (i == j) {
      return true;
    } else if (i > j) {
      return stored_in_lower(i, j);
    } else {
      assert(i < j);
      return stored_in_upper(i, j);
    }
  }

  // Get storage location {{{

  // The ~get_storage_location~ function tells us whether the (i, j)-th
  // coefficient is stored 'd': on the diagonal, ~'u'~: in the upper
  // triangle, ~'l'~: in the lower triangle, or ~'n'~: nowhere.


  char storage_location(int i, int j) const
  {
    if (i == j){
      return 'd';
    } else if (! stored(i, j)){
      return 'n';
    } else if (i < j) {
      return 'u';
    } else {
      assert (j < i);
      return 'l';
    }
  }

  // }}}

  // }}}

  struct not_found_error{}; 

  int get_position_in_upper(int i, int j) const
  {
    for (int m = 0; m < num_off_diag_elts / 2; m++) {
      if (owner[m] == i && neighbour[m] == j) {
        return m;
      }
    }
    throw not_found_error();
  }

  int get_position_in_lower(int i, int j) const
  {
    // Since the sparseness pattern is symmetric:
    return get_position_in_upper(j, i);
  }

  int index(int i, int j) const // {{{
  {

    /* {{{ Note we index everything from 0.  Note that, even if $A_{i,j}$ is
       0, it may still be stored, for example if i == j, or if $A_{j, i}$ is
       non-zero, or if $A_{i,j}$ used to be non-zero before someone modified A.

       But if $A_{i,j}$ is not stored then $A_{i,j}$ is definitely 0. }}} */

    char location = storage_location(i, j);

    if (location == 'n') {
      // not stored
      return -1000;
    } else if (location == 'd') {
      // diagonal
      return i;
    } else if (location == 'u') {
      // upper triangle
      return get_position_in_upper(i, j);
    } else {
      // lower triangle
      assert(location == 'l');
      return get_position_in_lower(i, j);
    }
  }
  // }}}

  void fill_in(int i, int j) // {{{
  {
    // We assume that (i, j) is not already in the sparseness pattern.
    assert(i != j);
    owner.push_back(i);
    neighbour.push_back(j);
    /* Update the */ matrix_band = std::max(matrix_band, abs(j - i));
    /* Update the */ num_off_diag_elts += 2; 
  }
  // Unlike with the compressed row format, we don't need to return the
  // position of the new entry, because here the position of the new
  // entry is just the end of the |vector|.
  // }}}

  void visualize_sparseness_pattern() const // {{{
  {
    std::vector<char> blank_row(n, '.');
    std::vector<std::vector<char>> grid(n, blank_row);

    for (int m = 0; m < num_off_diag_elts / 2; m++) {
      grid[owner[m]][neighbour[m]] = '#'; // A_{ij}
      grid[neighbour[m]][owner[m]] = '#'; // A_{ji}
    }

    for (int i = 0; i < n; i++) {
      grid[i][i] = '#';
    }

    // Print the sparseness pattern:
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	std::cout << grid[i][j];
      }
      std::cout << '\n';
    }
  }
  // }}}

  // Getters for the private data.{{{

  const std::vector<int> & get_owner() const { return owner; };
  const std::vector<int> & get_neighbour() const { return neighbour; }
  int get_n() const { return n; }
  int get_matrix_band () const { return matrix_band; }
  int get_num_off_diag_elts () const { return num_off_diag_elts; }
  // }}}
};

#endif // arrow_addresser_dot_H
