// -*- origami-fold-style: triple-braces; -*-

#ifndef compressed_row_addresser_dot_H
#define compressed_row_addresser_dot_H

// Compile with -D_GLIBCXX_DEBUG to get bounds checking

#include "unstructured-mesh.H"
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <list>
#include <cassert>

class compressed_row_addresser {
  // We assume that every row stores its diagonal coefficient.
  private:
    std::vector<int> row;
    std::vector<int> col;
    int band;
  public:

    // Constructor {{{

    // We do not assume that the cells or the neighbours are ordered.

    compressed_row_addresser(const unstructured_mesh &m)
    {
      row.resize(1 + m.cells.size());

      band = 0;

      int beginning_of_prev_row = 0;
      int size_of_prev_row = 0;

      for (cell c : m.cells) {
        int i = c.index;

        // Set |row[i]| to the position of the beginning of row i. {{{
        row[i] = beginning_of_prev_row + size_of_prev_row;
        // }}}

        // To |col|, add the indices of c's neighbours and c itself, in order. {{{

        std::list<int> s = c.neighbours;
        s.push_back(c.index);

        s.sort();

        for (int j : s) {
          col.push_back(j);
        }
        // }}}

        // Update the matrix band. {{{
        for (int j : c.neighbours) {
          band = std::max(band, abs(j - i));
        }
        // }}}

        beginning_of_prev_row = row[i];
        size_of_prev_row = 1 + c.neighbours.size();
      }

      // Add a ghost row. {{{
      row[m.cells.size()] = col.size();
      // A check {{{
      assert(col.size() == (size_t) (beginning_of_prev_row + size_of_prev_row));
      // }}}
      // }}}
    }
    // }}}

    // Resetting the floating point values is done by the matrix class.

    // Fill-in {{{

    // Assume that A_{i, j} is not yet stored in the matrix.

    int fill_in(int i, int j)
    {
      // Find the position in |col| where we should insert `j'. {{{
      int m;
      for (m = row[i]; m < row[i + 1]; m++) {
        assert(j != col[m]);
        if (j < col[m]) {
          break;
        }
      }
      // }}}

      // Insert j into col. {{{
      std::vector<int>::iterator iter = col.begin() + m;
      col.insert(iter, j);
      // }}}

      // Increment row[r] for r > i {{{
      for (int r = i+1; r < row.size(); r++){
        row[r]++;
      }
      // }}}

      // Update band {{{
      band = std::max(band, abs(j - i));
      // }}}

      // Return the position so that the matrix class can insert A_{i, j} into |coef| quickly.
      return m;
    }

    // The matrix class needs to insert A_{i, j} into |coef| by itself.
    // }}}

    int N () const { return row.size() - 1; } // Number of equations

    int num_off_diag_elts() const { return col.size() - N(); }

    // Get the index in |coef| from (i, j) {{{

    std::pair<bool, int> get_stored_and_index(int i, int j) const
    {
      // Search through row i for a coefficient with column value j.
      // If one is found, return |true| and its position within |row|.
      // Otherwise return |false| and some garbage integer.
      // Since col[i : i+1] is sorted, we exit early if j < col[m]

      for (int m = row[i]; m < row[i + 1]; m++) {
        if (j == col[m])
          return std::pair(true, m);
        else if (j < col[m])
          return std::pair(false, -2000);
      }
      return std::pair(false, -1000);
    }

    bool stored(int i, int j) const { return get_stored_and_index(i, j).first;}
    int index(int i, int j) const { return get_stored_and_index(i, j).second; }
    // }}}

    // Getters for the private data {{{
    const std::vector<int> & get_row() const { return row; }
    const std::vector<int> & get_col() const { return col; }
    int get_band() const { return band; }
    // }}}

    // Visualize sparseness pattern {{{
    void visualize_sparseness_pattern()
    {
      for (int i = 0; i < N(); i++) {
        // Visualize row i {{{
        std::vector<char> buffer(N(), ' ');
        // Fill buffer {{{
        for (int m = row[i]; m < row[i + 1]; m++) {
          int j = col[m];
          buffer[j] = '#';
        }
        // }}}
        // Print buffer {{{
        for (char c : buffer) {
          std::cout << c;
        }
        // }}}
        // }}}
        std::cout << std::endl;
      }
    }
    // }}}
};

#endif // compressed_row_addresser_dot_H
