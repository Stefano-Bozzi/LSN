CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio1.3 : esercizio1.3.o random.o
	$(CC) random.o esercizio1.3.o -o esercizio1.3	
esercizio1.3.o : esercizio1.3.cpp
	$(CC) -c esercizio1.3.cpp -o esercizio1.3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio1.3 seed.out
cleantxt :
	rm *.txt
