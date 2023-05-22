
CFLAGS = -c -Wall -g
LFLAGS = -g -o lodestar

all: lodestar clean

lodestar: lodestar.o ArgumentParser.o Engine.o
	g++ lodestar.o Engine.o $(LFLAGS)

lodestar.o: lodestar.cpp Engine.hpp
	g++ $(CFLAGS) lodestar.cpp

ArgumentParser.o: ArgumentPaser.cpp ArgumentParser.hpp
	g++ $(CFLAGS) ArgumentParser.cpp

Engine.o: Engine.cpp MatrixOperations.hpp MultidimensionalScaling.hpp Procrustes.hpp
	g++ $(CFLAGS) Engine.cpp

Procrustes.o: Procrustes.cpp Procrustes.hpp MatrixOperations.hpp NumericalRecipesInC.hpp
	g++ $(CFLAGS) Procrustes.cpp

MultidimensionalScaling.o: MultidimensionalScaling.cpp MultidimensionalScaling.hpp NumericalRecipesInC.hpp MatrixOperations.hpp
	g++ $(CFLAGS) MultidimensionalScaling.cpp

NumericalRecipesInC.o: NumericalRecipesInC.cpp NumericalRecipesInC.hpp
	g++ $(CFLAGS) NumericalRecipesInC.cpp

MatrixOperations.o: MatrixOperations.cpp MatrixOperations.hpp
	g++ $(CFLAGS) MatrixOperations.cpp

clean:
	rm -f *.o *~ *#