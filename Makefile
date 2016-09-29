
divergence: divergence.cpp divergence.hpp divergence.F90.o timer.o Makefile
	icc -O3 -restrict -vec-report=7 -xhost -Wall -std=c++11 divergence.cpp divergence.F90.o timer.o -o divergence -lrt

divergence.F90.o: divergence.F90 Makefile
	ifort -O3 -vec-report=7 -xhost -c -o divergence.F90.o divergence.F90

timer.o: timer/timer.cpp timer/timer.hpp
	icc -O3 -g -Wall -std=c++14 -c -o timer.o timer/timer.cpp

clean:
	rm -f divergence divergence.F90.o
