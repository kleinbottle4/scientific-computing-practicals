// -*- origami-fold-style: triple-braces; -*-

// Include files {{{
#define _GLIBCXX_DEBUG
#include <cmath>
#include <iostream>
#include <vector>
#include "mpi.h"
#include <cassert>
// }}}

bool prime(int N) // {{{
{
  // (To determine whether a number is prime, do not try anything
  // fancy; just check all possible odd divisors less than √N)
  if (N == 0) 
    return false;
  else if (N == 1) 
    return false;
  else if (N == 2) 
    return true;
  else if (N % 2 == 0) 
    return false;
  else 
    for (int d = 3; d <= sqrt(N); d += 2) 
      if (N % d == 0) 
	return false;
  return true;
  
}
// }}}

int ceiling(int a, int b) // {{{
{
  if (a % b == 0) {
    return a / b;
  } else {
    return (a/b) + 1;
  }
}
// }}}

#define DEBUG(x) std::cout << "debug rank " << rank << " " \
  << #x << " " << x << std::endl

// Constants
enum { DUMMY, TIME_TO_STOP, FIND_PRIMES };

int main()
{
  const int maximum = 10;
  const int block_size = 2;

  // Set up MPI {{{
  MPI_Init(NULL, NULL);
  int rank;
  int nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  // }}}

  // Some checks {{{
  assert(block_size <= maximum);
  assert(block_size >= 2);
  assert(nproc >= 2);
  // }}}

  if (rank == 0) {
    // Master code {{{
    int start = 0;
    int end = block_size;
    const int numblocks = ceiling(maximum, block_size);

    std::vector<int> primes;

    int number_of_stop_requests = 0; // Used as a check later.
    // Send the initial block limits to the other processes. {{{
    for (int r = 1; r <= nproc-1 && end < maximum; r++) {
      int startend[2] = {start, end};
      MPI_Send(startend, 2, MPI_INT, r, FIND_PRIMES, MPI_COMM_WORLD);
      start += block_size;
      end = std::min(end + block_size, maximum);
    }
    // }}}

    for (int numrecvd = 0; numrecvd < numblocks; numrecvd++) {
      DEBUG(numrecvd << start << end << number_of_stop_requests);
      // DO NOT assert(number_of_stop_requests < nproc - 1)

      // Receive a block of primes from any process, r. {{{
      // Count the number of new primes to add. {{{
      int count;
      MPI_Status status;
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_INT, &count);
      int r = status.MPI_SOURCE;
      // }}}
      // Use the tag to mark the start of the range.
      int tag = status.MPI_TAG;
      int oldsize = primes.size();
      primes.resize(oldsize + count);
      // Receive the new primes into the tail end of |primes|. {{{
      MPI_Recv((count == 0 ? NULL : &primes[oldsize]),
	       count, MPI_INT, r, tag, MPI_COMM_WORLD, &status);
      // }}}
      // }}}

      // Print the results. {{{
      std::cout << "Primes between" << tag << " and " << tag + block_size << ":\n";
      for (int j = oldsize; j < oldsize + count; j++) {
	std::cout << primes[j] << ", " << std::endl;
      }
      std::cout << std::endl;
      // }}}

      if (end < maximum) {
	// Increase |start| and |end| {{{
	start += block_size;
	end += block_size;
	if (end > maximum){ 
	  end = maximum;
	}
	// }}}
	// Send the limits of next block of naturals to that same process {{{
	int message[2] = {start, end};
	MPI_Send(message, 2, MPI_INT, r, FIND_PRIMES, MPI_COMM_WORLD);
	// }}}
      } else {
	assert(end == maximum);
	// Tell this process to stop. {{{
	MPI_Send(NULL, 0, MPI_INT, r, TIME_TO_STOP, MPI_COMM_WORLD);
	number_of_stop_requests++;
	// }}}
      }
    }
    assert(number_of_stop_requests == nproc - 1);
    // }}}
  }

  if (rank != 0) {
    // Slave code {{{
    int tag = DUMMY;
    while (tag != TIME_TO_STOP) {
      int startend[2];
      // Receive a message from the master. {{{
      MPI_Status status;
      MPI_Recv(startend, 2, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      tag = status.MPI_TAG;
      // }}}
      if (tag == FIND_PRIMES) {
	// Find primes in the interval requested by the master. {{{
	int start = startend[0];
	int end = startend[1];
	std::vector<int> primes;
	for (int j = start; j < end; j++) {
	  if (prime(j)) {
	    primes.push_back(j);
	  }
	}
	// }}}
	// Send the results to the master. {{{
	// The tag should be |start|.
	MPI_Send(primes.data(), primes.size(), MPI_INT, 0, start, MPI_COMM_WORLD);
	// }}}
      } else {
	assert(tag == TIME_TO_STOP);
      }
    }
    // }}}
  }

  MPI_Finalize();
  return 0;
}
