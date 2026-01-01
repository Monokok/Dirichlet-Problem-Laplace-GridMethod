#include "seq.cpp"
#include "mpi.cpp"


void main() {
	sequentialSolution();
	//parallelSolution(); //<- вместо вызова лучше скомпилировать сделать rename parallelSolution в main.cpp, собрать и запускать из-под консоли (стандарт MPI)
}