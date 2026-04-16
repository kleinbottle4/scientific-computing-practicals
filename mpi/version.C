#include <iostream>
#include "mpi.h"

int main() {
	int major, minor;
	MPI_Get_version(&major, &minor);
	std::cout << "Version is " << major << "." << minor << std::endl;
	return 0;
}
