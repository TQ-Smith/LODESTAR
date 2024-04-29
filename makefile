
# File: Makefile
# Date: 
# Author: TQ Smith
# Purpose: 

CC = gcc
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/lodestar: src/main.o
	$(CC) $(LFLAGS) bin/lodestar src/*.o lib/lapack/*.o -lz -lm -lpthread -lgfortran

src/main.o: src/SlidingWindow.o src/ProcrustesAnalysis.o
	$(CC) $(CFLAGS) src/main.c -o src/main.o

src/ProcrustesAnalysis.o: src/RealSymEigen.o
	$(CC) $(CFLAGS) src/ProcrustesAnalysis.c -o src/ProcrustesAnalysis.o

src/SlidingWindow.o: src/HaplotypeEncoder.o src/Window.o src/AlleleSharingDistance.o src/MultidimensionalScaling.o
	$(CC) $(CFLAGS) src/SlidingWindow.c -o src/SlidingWindow.o

src/MultidimensionalScaling.o: src/RealSymEigen.o
	$(CC) $(CFLAGS) src/MultidimensionalScaling.c -o src/MultidimensionalScaling.o

src/RealSymEigen.o: lib/lapack 
	$(CC) $(CFLAGS) src/RealSymEigen.c -o src/RealSymEigen.o

src/AlleleSharingDistance.o:
	$(CC) $(CFLAGS) src/AlleleSharingDistance.c -o src/AlleleSharingDistance.o

src/Window.o:
	$(CC) $(CFLAGS) src/Window.c -o src/Window.o

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o: src/RegionFilter.o
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

src/RegionFilter.o:
	$(CC) $(CFLAGS) src/RegionFilter.c -o src/RegionFilter.o

lib/lapack:
	$(CC) $(CFLAGS) lib/lapack/*.f -o lib/lapack

.PHONY: clean
clean:
	rm src/*.o bin/* lapack/*.o