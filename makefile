
CFLAGS = -c -Wall -g
LFLAGS = -g -o

bin/VCFGenotypeParser: src/VCFGenotypeParser.o
	gcc $(LFLAGS) bin/VCFGenotypeParser src/VCFGenotypeParser.o klib/kstring.o -lz

src/VCFGenotypeParser.o: klib/kstring.o
	gcc $(CFLAGS) src/VCFGenotypeParser.c -o src/VCFGenotypeParser.o

klib/kstring.o:
	gcc $(CFLAGS) klib/kstring.c -o klib/kstring.o

.PHONY: clean
clean:
	rm klib/*.o src/*.o bin/*