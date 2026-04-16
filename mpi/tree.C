#include "mpi.h"
#include <iostream>

int main()
{
  // Set up MPI {{{
  MPI_Init(NULL, NULL);
  int rank, n;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n);
  // }}}

  int v;
  
  // On process 0, initialize an integer v with value 0
  if (rank == 0) {
    v = 0;
  }

  // Send v from process 0 to process 1, and on process 1, add
  // 1. Display the value of v.

  // ...

  // Send v from process n to process n + 1, and on process n + 1, add
  // n. Display the value of v.

  for (int j = 0; j <= n-2; j++){
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == j){
      MPI_Send(&v, 1, MPI_INT, j+1, 100, MPI_COMM_WORLD);
    }
    if (rank == j+1){
      MPI_Recv(&v, 1, MPI_INT, j, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      v += rank;
      std::cout << "Process " << rank << " says v = " << v << ".\n";
    }
  }

  // Finally send this integer back to process zero and print out its
  // final value.

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == n-1){
    MPI_Send(&v, 1, MPI_INT, 0, 200, MPI_COMM_WORLD);
  }
  if (rank == 0){
    MPI_Recv(&v, 1, MPI_INT, n-1, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::cout << "I am process 0. v = " << v << ".\n";
  }

  // Clean up
  MPI_Finalize();

  return 0;
}
