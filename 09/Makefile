CC = g++
CFLAGS = -Wall -O3 --std=c++11

GA.exe : GA.o random.o lib.o
	$(CC) random.o lib.o GA.o -o GA.exe
GA.o : GA.cpp
	$(CC) -c GA.cpp -o GA.o $(CFLAGS)
lib.o : lib.cpp lib.h random.h
	$(CC) -c lib.cpp -o lib.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o GA.exe
cleandat :
	rm av_L2_Square.dat best_L2_Square.dat Square_final_Path_0.dat Square_first_Path_0.dat
