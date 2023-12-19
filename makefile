
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/HaplotypeTree: src/HaplotypeTree.o
	gcc $(LFLAGS) bin/HaplotypeTree src/HaplotypeTree.o src/VCFGenotypeParser.o klib/kstring.o -lz

src/HaplotypeTree.o: src/VCFGenotypeParser.o
	gcc $(CFLAGS) src/HaplotypeTree.c -o src/HaplotypeTree.o

src/VCFGenotypeParser.o: klib/kstring.o
	gcc $(CFLAGS) src/VCFGenotypeParser.c -o src/VCFGenotypeParser.o

klib/kstring.o:
	gcc $(CFLAGS) klib/kstring.c -o klib/kstring.o

.PHONY: clean
clean:
	rm klib/*.o src/*.o bin/*