
# File: Makefile
# Date: 9 May 2024
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Build LODESTAR. 

CC?=gcc
CFLAGS = -c -Wall -g -I lib
LFLAGS = -g -o

bin/lodestar: src/Main.o
	mkdir -p bin
	$(CC) $(LFLAGS) bin/lodestar src/*.o lib/lapack/*.o -lz -lm -lpthread -lgfortran

src/Main.o: src/SlidingWindow.o src/ProcrustesAnalysis.o src/Interface.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/ProcrustesAnalysis.o: src/Logger.o src/Matrix.o src/RealSymEigen.o
	$(CC) $(CFLAGS) src/ProcrustesAnalysis.c -o src/ProcrustesAnalysis.o

src/SlidingWindow.o: src/HaplotypeEncoder.o src/Window.o src/AlleleSharingDistance.o src/MultidimensionalScaling.o
	$(CC) $(CFLAGS) src/SlidingWindow.c -o src/SlidingWindow.o

src/Interface.o:
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

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

src/VCFLocusParser.o:
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

src/Matrix.o:
	$(CC) $(CFLAGS) src/Matrix.c -o src/Matrix.o

src/Logger.o:
	$(CC) $(CFLAGS) src/Logger.c -o src/Logger.o

.PHONY: lib/lapack
lib/lapack:
	gfortran $(CFLAGS) lib/lapack/*.f
	mv *.o lib/lapack

.PHONY: clean
clean:
	rm src/*.o bin/* lib/lapack/*.o