
CFLAGS = -c -Wall -g
LFLAGS = -g -o lodestar

all: lodestar

loadstar: lodestar.o Procrustes.o MultidimensionalScaling.o
	g++ lodestar.o Procrustes.o MultidimensionalScaling.o $(LFLAGS)

loadstar.o: lodestar.cpp Procrustes.h MultidimensionalScaling.h
	g++ $(CFLAGS) lodestar.cpp

Procrustes.o: Procrustes.cpp Procrustes.h MultidimensionalScaling.h NumericalRecipesInC.h
	g++ $(CFLAGS) Procrustes.cpp

MultidimensionalScaling.o: MultidimensionalScaling.cpp MultidimensionalScaling.h NumericalRecipesInC.h MatrixOperations.h
	g++ $(CFLAGS) MultidimensionalScaling.cpp

NumericalRecipesInC.o: NumericalRecipesInC.cpp NumericalRecipesInC.h
	g++ $(CFLAGS) NumericalRecipesInC.cpp

MatrixOperations.o: MatrixOperations.cpp MatrixOperations.h
	g++ $(CFLAGS) MatrixOperations.cpp

clean:
	rm -f lodestar *.o *~ *#