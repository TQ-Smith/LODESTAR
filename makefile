
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

src/Main.o: src/LODESTAR.o src/Interface.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interface.o: lib/kstring.o
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

src/LODESTAR.o: src/HaplotypeEncoder.o src/BlockList.o src/MatrixOperations.o
	$(CC) $(CFLAGS) src/LODESTAR.c -o src/LODESTAR.o

src/BlockList.o:
	$(CC) $(CFLAGS) src/BlockList.c -o src/BlockList.o

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	$(CC) $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o:
	$(CC) $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

src/MatrixOperations.o: lib/lapack 
	$(CC) $(CFLAGS) src/MatrixOperations.c -o src/MatrixOperations.o

lib/kstring.o:
	$(CC) $(CFLAGS) lib/kstring.c -o src/kstring.o

.PHONY: lib/lapack
lib/lapack:
	$(CC) $(CFLAGS) lib/lapack/*.f
	mv *.o lib/lapack

.PHONY: clean
clean:
	rm bin/* src/*.o lib/lapack/*.o lib/*.o