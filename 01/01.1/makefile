CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio1.1 : esercizio1.1.o random.o
	$(CC) random.o esercizio1.1.o -o esercizio1.1	
esercizio1.1.o : esercizio1.1.cpp
	$(CC) -c esercizio1.1.cpp -o esercizio1.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio1.1 seed.out
