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
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	int root = nproc - 1;

	double process_total = 0.0;
  
	int q = STEPS/nproc;
	int start = rank * q;
	int end   = start + q - 1;

	if (rank == nproc-1)
		end = STEPS-1;

	for(int i = start; i <= end; i++){
		double x = 10.0*((i+0.5)/(STEPS+1));
		process_total += exp( -pow(x, 2));
	}

        if (rank != 0) {
		MPI_Send(&process_total, 1, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
        }

	if (rank == 0){
		double grand_total = process_total;
		double buf;
		for (int i = 1; i <= nproc - 1; i++) {
			MPI_Recv(&buf, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			grand_total += buf;
		}

		std::cout << "Grand total =" << grand_total << ", integral = " << 10.0 * grand_total / STEPS << std::endl;

	}

	MPI_Finalize();

	return 0;
}
