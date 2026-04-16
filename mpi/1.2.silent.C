// -*- origami-fold-style: triple-braces; -*-

// Still broken!! 

// We want to re-write the Mandelbrot code so that it divides the grid
// into horizontal rows and divides them fairly between the processes,
// rather than just making `nproc` simple contiguous blocks.


#define _GLIBCXX_DEBUG

// Include files {{{
#include <cassert>
#include <iostream>
#include <fstream>
#include "mpi.h"
#include <vector>
#include <complex>
// }}}

// Grid class {{{


template <typename T>
class grid {

  // A two-dimensionally indexed vector with contiguous storage.  The
  // dimensions are set once, during construction. The only accessor is
  // operator()(int i, int j), which returns a reference to the i,j th
  // element of the grid. To get a pointer to the beginning of the i-th
  // row, use `&mygrid(i, 0)`.

private:
  const int number_of_rows;
  const int length_of_row;
  std::vector<T> v;
public:

  grid(int num, int len) : number_of_rows(num), length_of_row(len), v(num*len) {}

  T& operator()(int i, int j)
  {
    if (i >= number_of_rows){
      std::cerr << "Warning! i is " << i << " but number_of_rows is " << number_of_rows << ".\n";
    }
    if (j >= length_of_row){
      std::cerr << "Warning! j is " << j << " but length_of_row is " << length_of_row << ".\n";
    }
    if (i < 0) {
      std::cerr << "Warning i = " << i << " is negative.\n";
    }
    if (j < 0) {
      std::cerr << "Warning j = " << j << " is negative.\n";
    }

    return v.at(i*length_of_row + j);
  }

};
    
// }}}

// Use Rutter's definition of Mandel(z) {{{

int mandel(std::complex<double> z0, int iters)
{
  /* Repeats 320 iterations of the recurrence relation to 
     find the smallest integer, n, such that |z_n| < 2 
     
     If none is found, returns i = iters. This signifies 
     that z0 is not in Mandelbrot set */
  int i;
  std::complex<double> z;

  z = z0;
  for(i = 1; i < iters; i++){
    z = z*z + z0;
    if (abs(z) > 2.0)
      {
	break;
      }
  }
  return i;
}

// }}}

// Use Rutter's colouring-in function. {{{

int getColor(int iter, int totaliters, int index)
{
  // Returns r, g, or b color value given iter (result from calling
  // function mandel) The color map is a linear interpolation of the
  // colors vector defined below
  static std::vector<std::vector<int>> colors = {{68, 1, 84}, {32, 164, 134}, {68, 1, 84}, {57, 86, 140}, {31, 150, 139}, {115, 208, 85}, {253, 231, 37}, {42, 120, 142}, {255, 255, 255}};
  
  //Gradient region
  int no_gradients = colors.size() - 1;
  double var = (double)no_gradients * (1.0 - (double)iter / (double)totaliters);
  int gr = (int) var;
  if (gr > no_gradients - 1) gr = no_gradients - 1;
  if (gr < 0) gr = 0;

  return colors[gr][index] + (int) round((colors[gr+1][index] - colors[gr][index]) * (var - (double)gr));
}

// }}}

int main()
{
  // Declare some variables which are common to all processes. {{{

  double xmin, xmax, ymin, ymax;
  int resolution, numiters;
  int rank, nproc;

  // }}}

  // Initialize MPI. {{{

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  // }}}
  
  // On process 0, get input for the parameters in the same way as Rutter. {{{

  if (rank == 0) {
    // std::cout << "Enter resolution" << std::endl;
    // std::cin >> resolution;
    // do {
    //   std::cout << "Enter xmin" << std::endl;
    //   std::cin >> xmin;
    //   std::cout << "Enter xmax" << std::endl;
    //   std::cin >> xmax;
    //   std::cout << "Enter ymin" << std::endl;
    //   std::cin >> ymin;
    //   std::cout << "Enter ymax" << std::endl;
    //   std::cin >> ymax;
    //   if (xmin >= xmax || ymin >= ymax)
    //	{
    //	  std::cout << "Invalid domain boundaries." << std::endl;
    //	}
    // } while(xmin >= xmax && ymin >= ymax);

    // std::cout << "Enter iterations" << std::endl;
    // std::cin >> numiters;

    xmin = -1.0;
    xmax = 1.0;
    ymin = -1.0;
    ymax = 1.0;
    resolution = 4;
    numiters = 50;
  }
  // }}}

  // Broadcast the parameters from process 0 to all the other processes. {{{

  MPI_Bcast(&xmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ymin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ymin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&resolution, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&numiters, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // }}}

  // Synchronize all the processes. {{{

  MPI_Barrier(MPI_COMM_WORLD);

  // }}}

  // Calculate the number of rows that this process deals with. {{{

  // We use the division algorithm from GRM to deal with the cases
  // where the number of processes does not divide the number of rows
  // in the grid.

  const int q = resolution / nproc;
  const int r = resolution - q*nproc;

  // The number of rows this process computes for:
  const int numrows = (rank < r ? q+1 : q);
  // }}}

  // Store the answers for all the rows given to this process in a 2D grid.

  grid<int> results(numrows, resolution);

  // Iterate through all the row indices in [0, resolution) such that
  // n == rank (mod nproc)

  int k=0; // <- The position in `results` where the results for row i is stored.
  for (int i=rank; i<resolution; i += nproc){

    // Compute Mandel(point) for each point in the row and store in |results| {{{ 

    for (int j=0; j<resolution; j++){
      std::complex<double> z;
      z.real(xmin + j * ((xmax-xmin) / resolution));
      z.imag((ymax - i * ((ymax-ymin) / resolution)));
      assert(k < numrows);
      results(k, j) = mandel(z, numiters);
    }
    // }}}

    k++;
  }

  // Allocate space to store results on proc 0 {{{

  // In process 0, we need to allocate space for a large buffer to
  // store all the results. For the other processes, we don't need to
  // allocate anything, but we do allocate just one integer, so that
  // &all_the_results(0, 0) does not throw an error.

  grid<int> all_the_results((rank ? 1 : resolution),
			    (rank ? 1 : resolution));
  // }}}

  // Synchronize all the processes. {{{

  MPI_Barrier(MPI_COMM_WORLD);
  // }}}

  // All the computation is finished now.

  // Gather the results for rows 0, 1, 2, ..., nproc-1.
  // Then for rows nproc, nproc+1, nproc+2, ..., nproc+(nproc-1).
  // Etc. {{{

  // Recall q is floor(resolution/nproc);
  for (int n=0; n < q; n++){
    MPI_Gather(&results(n, 0),
	       resolution,
	       MPI_INT,
	       &all_the_results(0, 0),
	       resolution,
	       MPI_INT,
	       0,
	       MPI_COMM_WORLD);
  }
  // }}}
			    
  // Send the last remaining rows {{{

  // The first r processes have q+1 rows, and the remaining ones have
  // q rows. So for the first r processes we have one row of results
  // yet to send. We send them now to process 0 using point-to-point
  // communication.

  if (rank < r) {
    MPI_Send(&results(q, 0),
	     resolution,
	     MPI_INT,
	     0,
	     100,
	     MPI_COMM_WORLD);
  }

  if (rank == 0) {
    for (int n = 0; n < r; n++) {
      MPI_Recv(&all_the_results(q * nproc + n, 0),
	       resolution,
	       MPI_INT,
	       n,
	       100,
	       MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
    }
  }
  // }}}

  // Now Rank 0 plots the results to "mandelbrot.pam". {{{

  if (rank == 0){

    // std::ofstream img;
    // img.open("mandelbrot.pam", std::ios::binary);
    // img << "P6\n" << resolution << " " << resolution << " 255\n";
    // for (int i=0; i<resolution; i++){
    //   // Print one line of results.
    //   grid<unsigned char> line(resolution, 3);
    //   for (int j = 0; j < resolution; j++) {
    // 	for (int index : {0, 1, 2}){
    // 	  line(j, index) = getColor(all_the_results(i, j), numiters, index);
    // 	}
    //   }
    //   
    //   img.write((char *)(&line(0, 0)), 3 * resolution);
    // }
    // img.close();

    for (int i=0; i < resolution; i++){
      for (int j=0; j<resolution; j++){
	std::cout << all_the_results(i, j) << '\t';
      }
      std::cout << '\n';
    }
  }
  // }}}

  // Clean up MPI and exit. {{{

  MPI_Finalize();
  return 0;
  // }}}
}
