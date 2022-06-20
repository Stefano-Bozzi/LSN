CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio2.2 : esercizio2.2.o random.o Posizione.o
	$(CC) random.o Posizione.o esercizio2.2.o -o esercizio2.2	
esercizio2.2.o : esercizio2.2.cpp
	$(CC) -c esercizio2.2.cpp -o esercizio2.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
Posizione.o : Posizione.cpp Posizione.h
	$(CC) -c Posizione.cpp -o Posizione.o $(CFLAGS)
clean :
	rm *.o esercizio2.2 seed.out
cleantxt :
	rm *.txt
