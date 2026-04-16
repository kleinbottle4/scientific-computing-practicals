#include <iostream>
#include <math.h>
#include <iomanip>
#include "mpi.h"

#define STEPS 27720000

int main(int argc, char *argv[])
{
	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int root = size - 1;

	double process_total = 0.0;
  
	int q = STEPS/size;
	int start = rank * q;
	int end   = start + q - 1;

	if (rank == size-1)
		end = STEPS-1;

	for(int i = start; i <= end; i++){
		double x = 10.0*((i+0.5)/(STEPS+1));
		process_total += exp( -pow(x, 2));
	}

	double grand_total;

        MPI_Allreduce(&process_total, &grand_total, 1, MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);

	if (rank == root){
		std::cout << std::fixed << std::setprecision(16);
		std::cout << "Grand total = " << grand_total << ", integral = " << 10.0 * grand_total / STEPS << std::endl;
	}

        MPI_Finalize();

	return 0;
}
