#include <complex>
#include <iostream>
#include <fstream>
#include <vector>
#include "mpi.h"

#define DEBUG(s) std::cout << "===Rank " << rank << " : " << #s << " " << s << std::endl
#define DEBUG(s)

std::vector<std::vector<int>> colors = {{68, 1, 84}, {32, 164, 134}, {68, 1, 84}, {57, 86, 140}, {31, 150, 139}, {115, 208, 85}, {253, 231, 37}, {42, 120, 142}, {255, 255, 255}};

int mandel(std::complex<double> z0, int iters);

int getColor(int iter, int totaliters, int index);

int main(){
  double xmin, xmax, ymin, ymax;
  int res, iters;
  std::vector<int> row;
  std::ofstream img;

  int rank, size;
  std::vector<unsigned char> line;

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (!rank) {
    std::cout << "Enter resolution" << std::endl;
    std::cin >> res;
    do {
      std::cout << "Enter xmin" << std::endl;
      std::cin >> xmin;
      std::cout << "Enter xmax" << std::endl;
      std::cin >> xmax;
      std::cout << "Enter ymin" << std::endl;
      std::cin >> ymin;
      std::cout << "Enter ymax" << std::endl;
      std::cin >> ymax;
      if (xmin >= xmax || ymin >= ymax) {
	std::cout << "Invalid domain boundaries." << std::endl;
      }
    } while(xmin >= xmax && ymin >= ymax);

    std::cout << "Enter iterations" << std::endl;
    std::cin >> iters;

    row.resize(res);
    line.resize(3*res);

    img.open("para-mandel.pam", std::ios::binary);

    // Open output file and write width and height of image
    img << "P6\n" << res << " " << res << " 255\n";

  }

  MPI_Bcast( &res, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xmin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&xmax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ymin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ymax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&iters,1, MPI_INT, 0, MPI_COMM_WORLD);

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
      row_segment.at(j - start) = mandel(z, iters);
      // DEBUG(z);
      // DEBUG(row_segment[j - start]);
    }

    MPI_Gatherv(row_segment.data(),
		row_segment.size(),
		MPI_INT,
		(rank == 0 ? row.data() : NULL),
		sendcounts.data(),
		displs.data(),
		MPI_INT,
		0,
		MPI_COMM_WORLD);


    if (rank == 0){
      // Creates colours
      for(int j = 0; j < res; j++) {
	for (int index = 0; index < 3; index++) {
	  line[3 * j + index] = getColor(row[j], iters, index);
	}
      }

      /* Writes row data to file. (3 * res is the total number 
	 of characters needed for colours in one row) */
      img.write((char *) line.data(), 3 * res);
    }
  }

  if (!rank)
    img.close();

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

int getColor(int iter, int totaliters, int index)
{
  // Returns r, g, or b color value given iter (result from calling function mandel)
  // The color map is a linear interpolation of the colors vector defined above
  
  //Gradient region
  int no_gradients = colors.size() - 1;
  double var = (double)no_gradients * (1.0 - (double)iter / (double)totaliters);
  int gr = (int) var;
  if (gr > no_gradients - 1) gr = no_gradients - 1;
  if (gr < 0) gr = 0;

  return colors[gr][index] + (int) round((colors[gr+1][index] - colors[gr][index]) * (var - (double)gr));
}
