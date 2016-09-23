
divergence: divergence.cpp divergence.F90.o timer.o
	mpic++ -O3 -fp-model strict -fopenmp -g -Wall -std=c++14 divergence.cpp divergence.F90.o timer.o -o divergence

divergence.F90.o: divergence.F90
	mpif90 -O3 -fp-model strict -fopenmp -g -c -o divergence.F90.o divergence.F90

timer.o: timer/timer.cpp timer/timer.hpp
	mpic++ -O3 -g -Wall -std=c++14 -c -o timer.o timer/timer.cpp

clean:
	rm -f divergence divergence.F90.o
