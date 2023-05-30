
CFLAGS = -c -Wall -g
LFLAGS = -g -o lodestar

all: lodestar clean

lodestar: lodestar.o CommandLineArgumentParser.o Engine.o GenotypeFileParser.o Procrustes.o MultidimensionalScaling.o NumericalRecipesInC.o MatrixOperations.o
	g++ lodestar.o CommandLineArgumentParser.o Engine.o GenotypeFileParser.o Procrustes.o MultidimensionalScaling.o NumericalRecipesInC.o MatrixOperations.o $(LFLAGS)

lodestar.o: lodestar.cpp CommandLineArgumentParser.hpp Engine.hpp
	g++ $(CFLAGS) lodestar.cpp

Engine.o: Engine.cpp GenotypeFileParser.hpp Procrustes.hpp MultidimensionalScaling.hpp MatrixOperations.hpp
	g++ $(CFLAGS) Engine.cpp

CommandLineArgumentParser.o: CommandLineArgumentParser.cpp CommandLineArgumentParser.hpp
	g++ $(CFLAGS) CommandLineArgumentParser.cpp

GenotypeFileParser.o: GenotypeFileParser.cpp GenotypeFileParser.hpp
	g++ $(CFLAGS) GenotypeFileParser.cpp

Procrustes.o: Procrustes.cpp Procrustes.hpp NumericalRecipesInC.hpp MatrixOperations.hpp 
	g++ $(CFLAGS) Procrustes.cpp

MultidimensionalScaling.o: MultidimensionalScaling.cpp MultidimensionalScaling.hpp NumericalRecipesInC.hpp MatrixOperations.hpp 
	g++ $(CFLAGS) MultidimensionalScaling.cpp

NumericalRecipesInC.o: NumericalRecipesInC.cpp NumericalRecipesInC.hpp
	g++ $(CFLAGS) NumericalRecipesInC.cpp

MatrixOperations.o: MatrixOperations.cpp MatrixOperations.hpp
	g++ $(CFLAGS) MatrixOperations.cpp

clean:
	rm -f *.o *~ *#