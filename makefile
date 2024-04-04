
# File: Makefile
# Date: 
# Author: TQ Smith
# Purpose: 

CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/test: src/test.o
	gcc $(LFLAGS) bin/test src/*.o -lz -lm -lpthread

src/test.o: src/SlidingWindow.o
	gcc $(CFLAGS) src/test.c -o src/test.o

src/SlidingWindow.o: src/HaplotypeEncoder.o src/Window.o src/AlleleSharingDistance.o
	gcc $(CFLAGS) src/SlidingWindow.c -o src/SlidingWindow.o

src/AlleleSharingDistance.o:
	gcc $(CFLAGS) src/AlleleSharingDistance.c -o src/AlleleSharingDistance.o

src/Window.o:
	gcc $(CFLAGS) src/Window.c -o src/Window.o

src/HaplotypeEncoder.o: src/VCFLocusParser.o
	gcc $(CFLAGS) src/HaplotypeEncoder.c -o src/HaplotypeEncoder.o

src/VCFLocusParser.o:
	gcc $(CFLAGS) src/VCFLocusParser.c -o src/VCFLocusParser.o

.PHONY: clean
clean:
	rm src/*.o bin/*