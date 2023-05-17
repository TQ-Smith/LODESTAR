
CFLAGS = -c -Wall -g
LFLAGS = -g -o lodestar

all: lodestar clean

lodestar: lodestar.o Engine.o
	g++ lodestar.o Engine.o $(LFLAGS)

lodestar.o: lodestar.cpp Engine.h
	g++ $(CFLAGS) lodestar.cpp

Engine.o: Engine.cpp MatrixOperations.h MultidimensionalScaling.h Procrustes.h
	g++ $(CFLAGS) Engine.cpp

Procrustes.o: Procrustes.cpp Procrustes.h MatrixOperations.h NumericalRecipesInC.h
	g++ $(CFLAGS) Procrustes.cpp

MultidimensionalScaling.o: MultidimensionalScaling.cpp MultidimensionalScaling.h NumericalRecipesInC.h MatrixOperations.h
	g++ $(CFLAGS) MultidimensionalScaling.cpp

NumericalRecipesInC.o: NumericalRecipesInC.cpp NumericalRecipesInC.h
	g++ $(CFLAGS) NumericalRecipesInC.cpp

MatrixOperations.o: MatrixOperations.cpp MatrixOperations.h
	g++ $(CFLAGS) MatrixOperations.cpp

clean:
	rm -f *.o *~ *#