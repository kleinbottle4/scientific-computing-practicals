#define _GLIBCXX_DEBUG
#include <complex>
#include <vector>
#include "mpi.h"
#include <cstdlib>

int mandel(std::complex<double> z0, int iters);

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  int iters = atoi(argv[1]);
  int res   = atoi(argv[2]);
  double xmax = atof(argv[3]);
  double xmin = atof(argv[4]);
  double ymax = atof(argv[5]);
  double ymin = atof(argv[6]);

  std::vector<int> row;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> displs(size);
  std::vector<int> sendcounts(size);

  for (int r = 0; r < size; r++) {
    int s = r * res / size;
    int e = (r + 1) * res / size;
    sendcounts[r] = e - s;
    displs[r] = s;
  }

  for(int i = 0; i < res; i++) {
    // Calculates mandel() for each row in the domain

    int start = rank * res / size;
    int end = (rank + 1) * res / size;
    std::vector<int> row_segment(end - start);

    for(int j = start; j < end; j++){
      std::complex<double> z;
      z.real(xmin + j * ((xmax-xmin) / res));
      z.imag((ymax - i * ((ymax-ymin) / res)));
      row_segment[j - start] = mandel(z, iters);
    }

    MPI_Gatherv(row_segment.data(), row_segment.size(), MPI_INT, (rank == 0 ? row.data() : NULL), sendcounts.data(), displs.data(), MPI_INT, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
  return 0; 
}

int mandel(std::complex<double> z0, int iters)
{
  //  std::cout  << z0 << std::endl;
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
