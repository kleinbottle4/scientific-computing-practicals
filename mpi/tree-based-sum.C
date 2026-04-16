// Local Variables:
// tab-width: 4
// End:

#include "mpi.h"
#include <iostream>

#define DEBUG(s) std::cout << "rank " << rank << ": " << s << std::endl

int mylog2(int a)   
{
	int i = 0;
	while (a > 0) {
		a /= 2;
		i++;
	}
	return i - 1;
}

int mypow2(int a) 
{
	int i = 1;
	while (a >= 1) {
		i*=2;
		a--;
	}
	return i;
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int nproc;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
     
	int x=(1+rank)*(1+rank);
     
	int p = mylog2(nproc);
	while (p >= 0) {
		int m = mypow2(p); // 2^p
		DEBUG("p, m = " << p << " " << m);
		if (rank >= 2 * m) {
			DEBUG("Doing nothing.");
		} else if (rank >= m) {
			MPI_Send(&x, 1, MPI_INT, rank - m, 0, MPI_COMM_WORLD);
			DEBUG("sent" << x << " to " << rank-m);
		} else if (rank <= m - 1 && rank + m <= nproc - 1) {
			int y;
			MPI_Status status;
			MPI_Recv(&y, 1, MPI_INT, rank + m, 0, MPI_COMM_WORLD, &status);
			DEBUG("received" << y << " from " << rank+m);
			x += y;
		}
		p--;
	}

	if (rank == 0) {
		std::cout << "Sum is " << x << std::endl;
	}

	MPI_Finalize();

	return 0;
}
