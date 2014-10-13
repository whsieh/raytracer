CC = g++
CFLAGS = -g -Wall -std=c++11
TARGET = as2

all: $(TARGET)

default: as2

as2: as2.cpp lodepng.o parsing.h util.h
	$(CC) $(CFLAGS) -o as2 as2.cpp -I ./eigen

lodepng.o: lodepng.cpp
	$(CC) $(CFLAGS) -c lodepng.cpp

clean:
	rm -f *.o *.h.gch as2
	rm -rf *.dSYM
