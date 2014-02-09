CC = g++
CFLAGS = -std=c++11 -O3 -Wall

all: overlap sa_lcp_test

overlap: io.o sais.o lcp.o overlapper.o myers.o main.o
	$(CC) $(CFLAGS) -lm -pthread io.o sais.o lcp.o overlapper.o myers.o main.o -o overlap

sa_lcp_test: io.o sais.o lcp.o sa_lcp_test.o
	$(CC) $(CFLAGS) -lm io.o sais.o lcp.o sa_lcp_test.o -o sa_lcp_test

io.o: io.cpp io.hpp
	$(CC) $(CFLAGS) -c io.cpp

sais.o: sais.cpp sais.hpp
	$(CC) $(CFLAGS) -c sais.cpp

lcp.o: lcp.cpp lcp.hpp
	$(CC) $(CFLAGS) -c lcp.cpp

overlapper.o : overlapper.cpp overlapper.hpp ./myers/myers.h
	$(CC) $(CFLAGS) -c overlapper.cpp

myers.o: ./myers/myers.cpp ./myers/myers.h
	$(CC) $(CFLAGS) -c ./myers/myers.cpp

main.o: main.cpp io.hpp sais.hpp lcp.hpp overlapper.hpp
	$(CC) $(CFLAGS) -c main.cpp

sa_lcp_test.o: sa_lcp_test.cpp io.hpp sais.hpp lcp.hpp
	$(CC) $(CFLAGS) -c sa_lcp_test.cpp

clean:
	-rm *.o overlap sa_lcp_test
