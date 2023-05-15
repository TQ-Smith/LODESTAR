
CFLAGS = -c -Wall -g
LFLAGS = -g -o MultidimensionalScaling

all: MultidimensionalScaling

MultidimensionalScaling: MultidimensionalScaling.o NumericalRecipesInC.o MatrixOperations.o
	g++ MultidimensionalScaling.o NumericalRecipesInC.o MatrixOperations.o $(LFLAGS)

MultidimensionalScaling.o: MultidimensionalScaling.cpp NumericalRecipesInC.h MatrixOperations.h
	g++ $(CFLAGS) MultidimensionalScaling.cpp

NumericalRecipesInC.o: NumericalRecipesInC.cpp NumericalRecipesInC.h
	g++ $(CFLAGS) NumericalRecipesInC.cpp

MatrixOperations.o: MatrixOperations.cpp MatrixOperations.h
	g++ $(CFLAGS) MatrixOperations.cpp

clean:
	rm -f MultidimensionalScaling *.o *~ *#