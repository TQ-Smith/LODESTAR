#
# File: Makefile
# Date: 1 August 2023
# Author: TQ Smith
# Principal Investigator: Dr. Zachary Szpiech
# Purpose: Make and run LODESTAR and all unit tests.
#

# Possibly add optimization flag later.
CFLAGS = -std=c++11 -c -Wall -g
LFLAGS = -g -o

# Make LODESTAR and the unit tests.
.PHONY: all
all: bin/lodestar bin/unit_tests

# Make the LODESTAR executable.
bin/lodestar: utils src src/lodestar.cpp utils/CommandLineArgumentParser.hpp
	g++ src/lodestar.cpp src/*.o utils/LinearAlgebra/*.o $(LFLAGS) bin/lodestar

# Make the utils.
.PHONY: utils
utils: utils/LinearAlgebra/MatrixOperations.o utils/LinearAlgebra/MultidimensionalScaling.o utils/LinearAlgebra/NumericalRecipesInC.o

# Make MatrixOperations.
utils/LinearAlgebra/MatrixOperations.o: utils/LinearAlgebra/MatrixOperations.cpp utils/LinearAlgebra/MatrixOperations.hpp
	g++ $(CFLAGS) utils/LinearAlgebra/MatrixOperations.cpp -o utils/LinearAlgebra/MatrixOperations.o

# Make MultidimensionalScaling.
utils/LinearAlgebra/MultidimensionalScaling.o: utils/LinearAlgebra/MultidimensionalScaling.cpp utils/LinearAlgebra/MultidimensionalScaling.hpp utils/LinearAlgebra/NumericalRecipesInC.hpp
	g++ $(CFLAGS) utils/LinearAlgebra/MultidimensionalScaling.cpp -o utils/LinearAlgebra/MultidimensionalScaling.o

# Make NumericalRecipesInC.
utils/LinearAlgebra/NumericalRecipesInC.o: utils/LinearAlgebra/NumericalRecipesInC.cpp utils/LinearAlgebra/NumericalRecipesInC.hpp
	g++ $(CFLAGS) utils/LinearAlgebra/NumericalRecipesInC.cpp -o utils/LinearAlgebra/NumericalRecipesInC.o

# Make the source files.
.PHONY: src
src: src/Procrustes.o

# Make Procrustes.
src/Procrustes.o: src/Procrustes.cpp src/Procrustes.hpp utils/LinearAlgebra/MatrixOperations.hpp
	g++ $(CFLAGS) src/Procrustes.cpp -o src/Procrustes.o

# Make the unit tests for LODESTAR.
#	Note including the utils unit tests.
.PHONY: unit_tests
bin/unit_tests: bin/unit_tests/unit_test_procrustes

# Create the Procrustes unit test.
bin/unit_tests/unit_test_procrustes: unit_tests/unit_test_procrustes.cpp src/Procrustes.hpp utils/LinearAlgebra/MatrixOperations.hpp
	g++ unit_tests/unit_test_procrustes.cpp src/Procrustes.o utils/LinearAlgebra/MatrixOperations.o $(LFLAGS) bin/unit_tests/unit_test_procrustes

# Clean up object files.
.PHONY: clean
clean:
	rm bin/* bin/unit_tests/* src/*.o utils/LinearAlgebra/*.o