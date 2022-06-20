CC = g++
MPI = mpicxx
CFLAGS = -Wall -O3 --std=c++11

a.out : GA.o random.o lib.o
	$(MPI) random.o lib.o GA.o -o a.out
GA.o : GA.cpp
	$(MPI) -c GA.cpp -o GA.o $(CFLAGS)
lib.o : lib.cpp lib.h random.h
	$(MPI) -c lib.cpp -o lib.o $(CFLAGS)
random.o : random.cpp random.h
	$(MPI) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o a.out
