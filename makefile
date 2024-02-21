
# File: Makefile
# Date: 18 Janurary 2023
# Author: TQ Smith
# Purpose: Compiles Sliding Window Algorithm.

CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/test: src/test.o
	gcc $(LFLAGS) bin/test src/*.o klib/*.o -lz -lpthread -lm

src/test.o: src/SlidingWindow.o
	gcc $(CFLAGS) src/test.c -o src/test.o

src/SlidingWindow.o: src/HaplotypeEncoder.o src/Window.o src/Matrix.o
	gcc $(CFLAGS) src/SlidingWindow.c -o src/SlidingWindow.o

src/HaplotypeEncoder.o: src/VCFGenotypeParser.o
	gcc $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFGenotypeParser.o: klib/kstring.o
	gcc $(CFLAGS) src/VCFGenotypeParser.c -o src/VCFGenotypeParser.o

src/Window.o:
	gcc $(CFLAGS) src/Window.c -o src/Window.o

src/Matrix.o:
	gcc $(CFLAGS) src/Matrix.c -o src/Matrix.o

klib/kstring.o:
	gcc $(CFLAGS) klib/kstring.c -o klib/kstring.o

.PHONY: clean
clean:
	rm klib/*.o src/*.o bin/*