
divergence: divergence.cpp divergence.F90.o
	mpic++ -g -Wall -std=c++14 -fopenmp divergence.cpp divergence.F90.o -o divergence

divergence.F90.o: divergence.F90
	mpif90 -g -fopenmp -c -o divergence.F90.o divergence.F90

clean:
	rm -f divergence divergence.F90.o
