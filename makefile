
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
	$(CC) $(LFLAGS) bin/lodestar src/*.o lib/lapack/*.o lib/gsl/*.o -lz -lm -lpthread -lgfortran

src/Main.o: src/LODESTAR.o src/Interface.o
	$(CC) $(CFLAGS) src/Main.c -o src/Main.o

src/Interface.o: lib/kstring.o
	$(CC) $(CFLAGS) src/Interface.c -o src/Interface.o

src/LODESTAR.o: lib/gsl src/HaplotypeEncoder.o src/BlockList.o src/MatrixOperations.o
	$(CC) $(CFLAGS) -DHAVE_INLINE src/LODESTAR.c -o src/LODESTAR.o

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

.PHONY: lib/gsl
lib/gsl:
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/error.c -o lib/gsl/error.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/message.c -o lib/gsl/message.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/stream.c -o lib/gsl/stream.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/default.c -o lib/gsl/default.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/rng.c -o lib/gsl/rng.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/mt.c -o lib/gsl/mt.o
	$(CC) $(CFLAGS) -DHAVE_INLINE lib/gsl/types.c -o lib/gsl/types.o


.PHONY: clean
clean:
	rm bin/* src/*.o lib/lapack/*.o lib/*.o lib/gsl/*.o