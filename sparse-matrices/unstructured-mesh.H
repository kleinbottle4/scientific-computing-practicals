// file: unstructured_mesh.H
// Used by: The addresser and matrix classes
// Requires: Nothing

// This file defines a simple unstructured mesh class to be used by
// the sparse addressing classes and sparse matrix classes.

#ifndef unstructured_mesh_dot_H
#define unstructured_mesh_dot_H

#include <iostream>
#include <list>
#include <vector>

struct cell {
  int index;
  // A neighbour is another cell with a face which is in contact with
  // a face of this cell.

  // We only store the indices of the neighbours, to avoid recursion.

  std::list<int> neighbours; 

  // A real cell class would have more information, such as the cell
  // volume and the geometries of the faces.
};


struct unstructured_mesh {
  std::vector<cell> cells;
  int number_of_cells () const
  {
    return cells.size();
  }
};
#endif // unstructured_mesh_dot_H

