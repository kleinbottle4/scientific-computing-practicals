#include <iostream>
#include <vector>
#include "mpi.h"

int main()
{
    MPI_Init(NULL, NULL);
    int N;
    MPI_Comm_size(MPI_COMM_WORLD, &N);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> matrix;
    if (rank == 0){
	matrix.resize(N*N);
	for (int i=0; i<N; i++){
	    for (int j=0; j<N; j++){
		matrix[N*i + j] = (double) i*j+i+j+5;
	    }
	}
    }

    std::vector<double> row(N);

    MPI_Scatter(matrix.data(), N, MPI_DOUBLE, row.data(), N, MPI_DOUBLE, 0,
		MPI_COMM_WORLD);

    double sum = 0;
    for (double x : row) {
	sum += x;
    }

    std::cout << "Rank " << rank << " has a sum of " << sum << ".\n";

    std::vector<double> buf;
    if (rank == 0) {
	buf.resize(N);
    }

    MPI_Gather(&sum, 1, MPI_DOUBLE, buf.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
	double total_sum = 0;
	for (double x : buf) {
	    total_sum += x;
	}
	std::cout << "Rank 0 :: total sum = " << total_sum << std::endl;
    }

    MPI_Finalize();

    return 0;
}
