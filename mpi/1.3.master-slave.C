// -*- origami-fold-style: triple-braces; -*-

#define _GLIBCXX_DEBUG // To get bounds checking on vectors

// Header files {{{
#include "mpi.h"
#include <vector>
#include <cassert>
#include <complex>
#include <iostream>
#include <fstream>
// }}}

enum { CALCULATE, FINISH, DUMMY, RESULTS };

// Simple grid class {{{
template <typename T>
class grid {
  // n rows and m columns.
  // Stored row by row like a 2D C-style array.
  int n, m;
  std::vector<T> v;
public:
  grid(int n, int m)
  {
    grid::n = n;
    grid::m = m;
    v.resize(n * m);
  }
  T &operator()(int i, int j)
  {
    assert(0 <= i && i < n && 0 <= j && j < m);
    return v[m * i + j];
  }
};
// }}}

int mandel(std::complex<double> z0, int iters) // {{{
{
  int i;
  std::complex<double> z;

  z = z0;
  for(i = 1; i < iters; i++){
    z = z*z + z0;
    if (abs(z) > 2.0)
      break;
  }
  return i;
}
// }}}

// Colours {{{

std::vector<std::vector<int>> colors = {
  {68, 1, 84},    {32, 164, 134}, {68, 1, 84},
  {57, 86, 140},  {31, 150, 139}, {115, 208, 85},
  {253, 231, 37}, {42, 120, 142}, {255, 255, 255}
};

int getColor(int iter, int totaliters, int index)
{
  // Returns r, g, or b color value given iter (result from calling
  // function mandel) The color map is a linear interpolation of the
  // colors vector defined above
  
  //Gradient region
  int no_gradients = colors.size() - 1;
  double var = (double)no_gradients * (1.0 - (double)iter / (double)totaliters);

  int gr = (int) var;
  if (gr > no_gradients - 1) gr = no_gradients - 1;
  if (gr < 0) gr = 0;

  return colors[gr][index] + (int)round((colors[gr + 1][index] -
    colors[gr][index]) * (var - (double)gr));

}

// }}}

int main(int argc, char **argv)
{
  // Set up MPI {{{
  MPI_Init(&argc, &argv);
  int rank, n;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  // }}}

  // Parameters
  double xmin, xmax, ymin, ymax;
  int res, iters;

  if (rank == 0) {
    // Get parameters from user {{{
    res  = 40; // atoi(argv[1]);
    xmin = -2.0; // atof(argv[2]);
    xmax = 2.0;  // atof(argv[3]);
    ymin = -2.0; // atof(argv[4]);
    ymax = 2.0; // atof(argv[5]);
    iters = 20; // atoi(argv[6]);
    if (xmin >= xmax || ymin >= ymax) {
      std::cout << "Invalid domain boundaries." << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    // }}}
  }

  // Broadcast the parameters {{{
  MPI_Bcast(&res ,  1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&xmin,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xmax,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ymin,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ymax,  1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&iters, 1, MPI_INT,    0, MPI_COMM_WORLD);
  // }}}

  if (rank == 0) {
    // Set up the grid {{{
    grid<int> g(res, res);
    // }}}

    // Send initial row indices {{{

    // If the number of processes is more than the number of rows,
    // then the main loops do not get executed.
    
    int numsent = std::min(n - 1, res);
    for (int i=0; i < numsent; i++) {
      int rowidx = i;
      MPI_Send(&rowidx, 1, MPI_INT, i+1, CALCULATE, MPI_COMM_WORLD);
    }
    // }}}
      
    // Main loop for master {{{
    int numrecvd = 0;

    std::vector<int> row(res);
    while (numrecvd < res) {
      // A check {{{
      assert(numsent >= numrecvd);
      // }}}

      // Receive a calculated row from any other process {{{
      MPI_Status st;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
      // The tag is the row index.
      int i = st.MPI_TAG;
      int recvrank = st.MPI_SOURCE;
      MPI_Recv(&g(i, 0), res, MPI_INT, recvrank, i, MPI_COMM_WORLD, &st);
      // A check {{{
      assert(st.MPI_TAG == i);
      assert(st.MPI_SOURCE == recvrank);
      // }}}

      // }}}
      numrecvd++;

      if (numsent < res){
	// Send a new index (numsent-1+1) to the same process. {{{
	MPI_Send(&numsent, 1, MPI_INT, recvrank, CALCULATE, MPI_COMM_WORLD);
	// }}}
	numsent++;
      } else {
	// A check {{{
	assert(numsent == res);
	// }}}
	// Tell the slave to finish {{{
	int dummy;
	MPI_Send(&dummy, 1, MPI_INT, recvrank, FINISH, MPI_COMM_WORLD);
	// }}}
      }
    }
    // }}}

    MPI_Barrier(MPI_COMM_WORLD);

    // Print the results {{{
    std::ofstream img;
    img.open("mandel.pam", std::ios::binary);
    img << "P6\n" << res << " " << res << " 255\n";
    for (int i=0; i < res; i++){
      grid<unsigned char> line(res, 3);
      // std::vector<unsigned char> line(3*res);
      for (int j = 0; j < res; j++){
	for (int index : {0, 1, 2}) {
	  line(j, index) = getColor(g(i, j), iters, index);
	  // line[3*j + index] = getColor(g(i, j), iters, index);
	}
      }
      img.write((char *) &line(0, 0), 3*res);
      // img.write((char *) line.data(), 3*res);
    }
    img.close();
    // }}}

    /*
    // Debug : print results {{{
    for (int i = 0; i < res; i++){
      for (int j = 0; j < res; j++){
	std::cout << (g(i, j) == iters ? '#' : '.');
      }
      std::cout << std::endl;
    }
    // }}}
    */

  } else {
    int numrowscomputed = 0;
    // Main loop for slaves {{{

    while (true) {
      int i;
      MPI_Status st;

      // Receive  a message from the master {{{
      MPI_Recv(&i, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
      // }}}

      if (st.MPI_TAG == FINISH) {
	break;
      } else {
	assert(st.MPI_TAG == CALCULATE);
	// i is the index of the row to calculate.
	// Calculate Mandel(row i) {{{
	std::vector<int> row(res);
	for (int j = 0; j < res; j++) {
	  std::complex<double> z(xmin + j * (xmax - xmin) / res,
				 ymin + i * (ymax - ymin) / res);
	  row[j] = mandel(z, iters);
	}
	// }}}
	numrowscomputed++;
	// Send the results back to the master with tag i. {{{
	MPI_Send(row.data(), res, MPI_INT, 0, i, MPI_COMM_WORLD);
	// }}}
      }
    }
    // }}}
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "Process " << rank << " computed "
	      << numrowscomputed << " rows.\n";
  }

  MPI_Finalize();
  return 0;
}
