
# File: Makefile
# Date: 9 May 2024
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build LODESTAR. ONLY compiles with gcc. NOT clang.

CC = gcc-13
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/lodestar: src/Main.o
	$(CC) $(LFLAGS) bin/lodestar src/*.o lib/lapack/*.o -lz -lm -lpthread -lgfortran

src/Main.o: src/Logger.o src/SlidingWindow.o src/ProcrustesAnalysis.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

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

src/VCFLocusParser.o: src/RegionSet.o
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

src/RegionSet.o:
	$(CC) $(CFLAGS) src/RegionSet.c -o src/RegionSet.o

src/Logger.o:
	$(CC) $(CFLAGS) src/Logger.c -o src/Logger.o

lib/lapack:
	$(CC) $(CFLAGS) lib/lapack/*.f -o lib/lapack

.PHONY: clean
clean:
	rm src/*.o bin/* lapack/*.o