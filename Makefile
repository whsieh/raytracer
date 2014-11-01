CC = g++
CFLAGS = -g -Wall -O2 -std=c++0x
TARGET = trace

all: $(TARGET)

default: trace

trace: trace.cpp lodepng.o parsing.h util.h
	$(CC) $(CFLAGS) -o trace trace.cpp -I ./lib/eigen

lodepng.o: lib/lodepng.cpp
	$(CC) $(CFLAGS) -c lib/lodepng.cpp

clean:
	rm -f *.o *.h.gch trace
	rm -rf *.dSYM
