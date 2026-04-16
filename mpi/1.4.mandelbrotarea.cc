#include <complex>
#include <iostream>
#include "mpi.h"

#define DEBUG(s) std::cout << "===Rank " << rank << " : " << #s << " " << s << std::endl

int mandel(std::complex<double> z0, int iters);
bool inmandelbrotset(std::complex<double> z, int iters);

int main(int argc, char ** argv)
{
  double xmin = atof(argv[1]);
  double xmax = atof(argv[2]);
  double ymin = atof(argv[3]);
  double ymax = atof(argv[4]);
  int res   = atoi(argv[5]);
  int iters = atoi(argv[6]);

  int rank, size;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int numberofpointsinside = 0;
  int start = (rank * res) / size;
  int end = ((rank + 1) * res) / size;
  for(int i = start; i < end; i++) {
    for(int j = 0; j < res; j++){
      std::complex<double> z;
      z.real(xmin + j * ((xmax-xmin) / res));
      z.imag((ymax - i * ((ymax-ymin) / res)));
      if (inmandelbrotset(z, iters)) {
	numberofpointsinside++;
      }
    }
  }

  const double numberofpoints = (end - start) * res;
  double myymin = ymin + start * ((ymax - ymin) / res);
  double myymax = ymin + end   * ((ymax - ymin) / res);
  const double areaofmyregion = (xmax-xmin)*(myymax-myymin);
  // Cannot make |partialareaofmandelbrotset| const due to MPI.
  double area_of_my_mandelbrot_subset = areaofmyregion *
    (numberofpointsinside / numberofpoints);

  DEBUG(areaofmyregion);
  DEBUG(area_of_my_mandelbrot_subset);
  DEBUG(numberofpoints);
  DEBUG(numberofpointsinside);

  double areaofwholemandelbrotset;

  MPI_Reduce(&area_of_my_mandelbrot_subset, &areaofwholemandelbrotset, 1,
	     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0){
    std::cout << "The area of the whole mandelbrot set is about "
	      << areaofwholemandelbrotset << ".\n";
  }

  MPI_Finalize();
  return 0; 
}

int mandel(std::complex<double> z0, int iters)
{
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

bool inmandelbrotset(std::complex<double> z, int iters)
{
  return mandel(z, iters) >= iters;
}
