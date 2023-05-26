
CFLAGS = -c -Wall -g
LFLAGS = -g -o lodestar

all: lodestar clean

lodestar: lodestar.o CommandLineArgumentParser.o Engine.o
	g++ lodestar.o Engine.o $(LFLAGS)

lodestar.o: lodestar.cpp Engine.hpp
	g++ $(CFLAGS) lodestar.cpp

CommandLineArgumentParser.o: CommandLineArgumentParser.cpp CommandLineArgumentParser.hpp
	g++ $(CFLAGS) CommandLineArgumentParser.cpp

Engine.o: Engine.cpp GenotypeFileParser.hpp MatrixOperations.hpp MultidimensionalScaling.hpp Procrustes.hpp
	g++ $(CFLAGS) Engine.cpp

GenotypeArgumentParser.o: GenotypeFileParser.cpp GenotypeFileParser.hpp
	g++ $(CFLAGS) GenotypeFileParser.cpp

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