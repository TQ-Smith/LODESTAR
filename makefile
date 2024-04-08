
# File: Makefile
# Date: 
# Author: TQ Smith
# Purpose: 

CC = gcc-13
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/test: src/test.o
	$(CC) $(LFLAGS) bin/test src/*.o lib/lapack/*.o -lz -lm -lpthread

src/test.o: src/SlidingWindow.o src/MultidimensionalScaling.o
	$(CC) $(CFLAGS) src/test.c -o src/test.o

src/SlidingWindow.o: src/HaplotypeEncoder.o src/Window.o src/AlleleSharingDistance.o
	$(CC) $(CFLAGS) src/SlidingWindow.c -o src/SlidingWindow.o

src/MultidimensionalScaling.o: lib/lapack
	$(CC) $(CFLAGS) src/MultidimensionalScaling.c -o src/MultidimensionalScaling.o

src/AlleleSharingDistance.o:
	$(CC) $(CFLAGS) src/AlleleSharingDistance.c -o src/AlleleSharingDistance.o

src/Window.o:
	$(CC) $(CFLAGS) src/Window.c -o src/Window.o

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o:
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

lib/lapack:
	gfortran $(CFLAGS) lib/lapack/*.f

.PHONY: clean
clean:
	rm src/*.o bin/* lapack/*.o