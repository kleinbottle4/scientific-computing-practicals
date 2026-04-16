#include <iostream>


int main()
{
    int rank;
    int size;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> matrix(N*N);
    for (int i=0; i<N; i++){
	for (int j=0; j<N; j++){
	    matrix[i][j] = i*j;
