// -*- origami-fold-style: triple-braces; -*-

/*

 */

// Header files {{{
#define _GLIBCXX_DEBUG
#include <fstream>
#include <math.h>
#include "mpi.h"
#include <vector>
#include <iostream>
#include <cassert>
// }}}

// Constants {{{
#define SIZE 400
enum { SENDING_BOTTOM_ROW, SENDING_TOP_ROW, SENDING_RESULTS };
// }}}

struct rectangle { // {{{
                   // 2D array class with variable starting index for the row.
  int start;
  int end;
  int height;
  int width;
  std::vector<double> v;

  rectangle(int start, int end, int width)
  {
    rectangle::start = start;
    rectangle::end = end;
    height = end - start;
    rectangle::width = width;
    v.resize(height * width);
  }

  double& operator()(int i, int j)
  {
    i -= start;
    if (! (i >= start)) {
      std::cout << i << " " << j << " " << start << ".\n";
    assert(i >= start);
    }
    assert(i < end);
    assert(j >= 0);
    assert(j < width);
    return v[(i - start)*width + j];
  }
};
// }}}

int main()
{
  // Set up MPI {{{
  MPI_Init(NULL, NULL);
  int nproc;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // }}}

  // Set up parameters {{{
  // i is the row.  j is the column. height is the height. niter is the
  // number of iterations. n is the current iteration number.

  int height=SIZE+2;
  int niter=4000;
  // }}}

  // Calculate the start and end positions for this process {{{
  int start = rank * height / nproc;
  int end = (rank + 1) * height / nproc;
  // }}}

  // Allocate the blocks of the grids which this process deals with. {{{
  rectangle g1(start - 1, end + 1, SIZE);
  rectangle g2(start, end, SIZE);
  // }}}

  // Fill the old grid with zeroes. {{{
  for (int i = start - 1; i < end + 1; i++) {
    for (int j = 0; j < SIZE; j++) {
      g1(i, j) = 0.0;
    }
  }
  // }}}

  // Apply the boundary conditions {{{
  if (rank == 0) {
    // Apply the top boundary condition to the top ghost row.
    for (int j = SIZE/8; j < 7*SIZE/8; j++) {
      g1(start - 1, j) = 1.0;
    }
  } else if (rank == nproc - 1) {
    // Apply the bottom boundary condition to the bottom ghost row.
    for (int j = 0; j < SIZE; j++) {
      g1(end, 0) = 0.0;
    }
  }
  // }}}

  // Set up the buffer for MPI_Bsend {{{
  int pack_size;
  MPI_Pack_size(SIZE, MPI_DOUBLE, MPI_COMM_WORLD, &pack_size);
  int num_sends = 2;
  int size_of_buffer = num_sends * (MPI_BSEND_OVERHEAD + pack_size); 
  double *buffer = (double *) malloc(size_of_buffer);
  MPI_Buffer_attach(buffer, size_of_buffer);
  // }}}

  for (int n = 0; n < niter; n++) {
    // Synchronize all the processes. {{{
    MPI_Barrier(MPI_COMM_WORLD);
    // }}}

    // Move information at boundaries downwards. {{{

    // Note that I use the convention that (0, 0) is at the top left
    // and that process 0 takes the top block.

    // Send my bottom row to the top ghost row of the process below me. {{{
    MPI_Bsend(&g1(end - 1, 0), SIZE, MPI_DOUBLE,
	     (rank == nproc - 1 ? MPI_PROC_NULL : rank + 1),
	     SENDING_BOTTOM_ROW, MPI_COMM_WORLD);
    // }}}
    // Get the row above my block and store it in the top ghost row. {{{
    MPI_Recv(&g1(start - 1, 0), SIZE, MPI_DOUBLE,
	     (rank == 0 ? MPI_PROC_NULL : rank - 1),
	     SENDING_BOTTOM_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // }}}
    // }}}

    // Move information at boundaries upwards. {{{

    // Send my top row to the bottom ghost row of the process above me. {{{
    MPI_Bsend(&g1(start, 0), SIZE, MPI_DOUBLE,
	     (rank == 0 ? MPI_PROC_NULL : rank - 1),
	     SENDING_TOP_ROW, MPI_COMM_WORLD);

    // }}}

    // Get the row below my block and store it in the bottom ghost row. {{{
    MPI_Recv(&g1(end, 0), SIZE, MPI_DOUBLE,
	     (rank == nproc-1 ? MPI_PROC_NULL : rank + 1),
	     SENDING_TOP_ROW, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // }}}
    // }}}

    // Apply the Laplacian numerical stencil to my block. {{{

    for (int i = start; i < end; i++) {
      // Left:
      g2(i, 0) = (1.0/3) * (g1(i, 1) + g1(i-1, 0) + g1(i+1, 0));
      // Right:
      g2(i, SIZE-1) = (1.0/3) * (g1(i, SIZE-2) + g1(i-1, SIZE-1) + g1(i+1, SIZE-1));
      // Middle:
      for (int j = 1; j < SIZE - 1; j++) {
        g2(i, j) = 0.25 * (g1(i-1, j) + g1(i+1, j) + g1(i, j-1) + g1(i, j+1));
      }
    }

    // }}}

    // Get rid of the old values. {{{
    for (int i = start; i < end; i++) {
      for (int j = 0; j < SIZE; j++) {
        g1(i, j) = g2(i, j);
      }
    }
    // Don't worry about the ghost values. They get overwritten in the
    // next iteration of the loop.

    // }}}
  }

  // Free the MPI_Bsend buffer {{{
  MPI_Buffer_detach(&buffer, &size_of_buffer);
  free(buffer);
  // }}}

  // Synchronize all the processes. {{{
  MPI_Barrier(MPI_COMM_WORLD);
  // }}}

  // Send the results back to process 0 in rank order {{{
  if (rank == 0) {
    rectangle biggrid(0, SIZE, SIZE);
    // Receive results into ~biggrid~ {{{
    for (int r = 0; r <= nproc-1; r++) {
      int s = r * SIZE / nproc;
      // int e = (r+1) * SIZE / nproc;
      MPI_Recv(&biggrid(s, 0), SIZE, MPI_DOUBLE, r,
	       SENDING_RESULTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // }}}

    // Print results to =laplace.pam= {{{
    std::ofstream img;
    img.open("laplace.pam");
    img << "P2\n" << SIZE << " " << SIZE << " 255\n";
    for (int i = 1; i < SIZE + 1; i++) {
      for (int j = 0; j < SIZE; j++) {
        img << (int) (255*biggrid(i, j)) << "\n";
      }
    }
    // }}}
  } else {
    // Send the results but not the ghost rows to rank 0. {{{
    MPI_Send(&g1(start, 0), (end - start) * SIZE, MPI_DOUBLE,
	     0, SENDING_RESULTS, MPI_COMM_WORLD);
    // }}}
  }
  // }}}

  // Clean up {{{
  MPI_Finalize();
  // }}}

  return 0;
}
