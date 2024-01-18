
# File: Makefile
# Date: 18 Janurary 2023
# Author: TQ Smith
# Purpose: Compiles Sliding Window Algorithm.

CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/SlidingWindow: src/SlidingWindow.o
	gcc $(LFLAGS) bin/SlidingWindow src/SlidingWindow.o src/HaplotypeEncoder.o src/VCFGenotypeParser.o klib/kstring.o -lz

src/SlidingWindow.o: src/HaplotypeEncoder.o
	gcc $(CFLAGS) src/SlidingWindow.c -o src/SlidingWindow.o

src/HaplotypeEncoder.o: src/VCFGenotypeParser.o
	gcc $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFGenotypeParser.o: klib/kstring.o
	gcc $(CFLAGS) src/VCFGenotypeParser.c -o src/VCFGenotypeParser.o

klib/kstring.o:
	gcc $(CFLAGS) klib/kstring.c -o klib/kstring.o

.PHONY: clean
clean:
	rm klib/*.o src/*.o bin/*