CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio3 : esercizio3.o random.o
	$(CC) random.o esercizio3.o -o esercizio3	
esercizio3.o : esercizio3.cpp
	$(CC) -c esercizio3.cpp -o esercizio3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio3 seed.out
cleantxt :
	rm *.txt
