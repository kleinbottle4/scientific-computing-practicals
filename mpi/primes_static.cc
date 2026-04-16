// -*- origami-fold-style: triple-braces; -*-

// Primes calculation

// Task {{{

// This question asks you to write a parallel program that computes
// the primes between 1 and 10000000, in two ways.  In both versions,
// the prime numbers, as well as the number of primes found by each
// process should be sent to rank 0, who will eventually print the
// total number of primes found.

// ** Version 1: Static workload

// First, divide the workload (approximately) equally between the
// processes.

// Hint: Use MPI Get count to determine the number of primes
// that each process has computed in each range.

// }}}

// Header files {{{
#define _GLIBCXX_DEBUG
#include <cmath>
#include <iostream>
#include <vector>
#include "mpi.h"
// }}}


bool prime(int N) // {{{
{
  // (To determine whether a number is prime, do not try anything
  // fancy; just check all possible odd divisors less than √N)
  if (N == 1) {
    return false;
  } else if (N == 2) {
    return true;
  } else if (N % 2 == 0) {
    return false;
  } else {
    for (int d = 3; d <= sqrt(N); d += 2) {
      if (N % d == 0) {
	return false;
      }
    }
    return true;
  }
}
// }}}


#define MAXIMUM 10000 // 10000

int main()
{
  // Set up MPI {{{
  MPI_Init(NULL, NULL);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // }}}

  // Calculate start and end points for this process. {{{
  int start = rank * (MAXIMUM / size);
  int end = (rank + 1) * (MAXIMUM / size);
  // We operate on [start, end).
  // }}}

  // Find all the primes in this process's range. {{{
  std::vector<int> primes;
  for (int n = start; n < end; n++) {
    if (prime(n)) {
      primes.push_back(n);
    }
  }
  // }}}

  if (rank != 0) {
    // Send the results back to process 0. {{{
    MPI_Send(primes.data(), primes.size(), MPI_INT, 0, 100, MPI_COMM_WORLD);
    // }}}
  }

  if (rank == 0) {
    // Receive primes from all processes, in order. {{{
    for (int r = 1; r <= size-1; r++) {
      int oldsize = primes.size();
      // In 'count', store number of primes from process r. {{{
      MPI_Status st;
      MPI_Probe(r, 100, MPI_COMM_WORLD, &st);
      int count;
      MPI_Get_count(&st, MPI_INT, &count);
      // }}}
      primes.resize(primes.size() + count);
      // Receive the primes from process r into the tail end of
      // allprimes. {{{
      MPI_Recv(&primes[oldsize], count, MPI_INT, r, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // }}}
    }
    // }}}
  }

  if (rank == 0) {
    // Print out the primes, in order. {{{
    for (int j = 0; j < primes.size(); j++) {
      std::cout << primes[j] << ", ";
      if (j % 20 == 0) {
	std::cout << "\n";
      }
    }
    // }}}
  }

  MPI_Finalize();
  return 0;

}

