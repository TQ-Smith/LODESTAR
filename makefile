
# File: Makefile
# Date: 
# Author: TQ Smith
# Purpose: 

CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/test: src/HaplotypeEncoder.o
	gcc $(LFLAGS) bin/test src/*.o -lz 

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	gcc $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o:
	gcc $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

.PHONY: clean
clean:
	rm src/*.o bin/*