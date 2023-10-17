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
utils/LinearAlgebra/MatrixOperations.o:
	g++ $(CFLAGS) utils/LinearAlgebra/MatrixOperations.cpp -o utils/LinearAlgebra/MatrixOperations.o

# Make MultidimensionalScaling.
utils/LinearAlgebra/MultidimensionalScaling.o:
	g++ $(CFLAGS) utils/LinearAlgebra/MultidimensionalScaling.cpp -o utils/LinearAlgebra/MultidimensionalScaling.o

# Make NumericalRecipesInC.
utils/LinearAlgebra/NumericalRecipesInC.o:
	g++ $(CFLAGS) utils/LinearAlgebra/NumericalRecipesInC.cpp -o utils/LinearAlgebra/NumericalRecipesInC.o

# Make the source files.
.PHONY: src
src: src/Procrustes.o src/VCFParser.o src/SlidingWindow.o

# Make Procrustes.
src/Procrustes.o: src/Procrustes.cpp src/Procrustes.hpp utils/LinearAlgebra/MatrixOperations.hpp
	g++ $(CFLAGS) src/Procrustes.cpp -o src/Procrustes.o

# Make VCFParser.
src/VCFParser.o: src/VCFParser.cpp src/VCFParser.hpp
	g++ $(CFLAGS) src/VCFParser.cpp -o src/VCFParser.o

# Make SlidingWindow.
src/SlidingWindow.o: src/SlidingWindow.cpp src/SlidingWindow.hpp
	g++ $(CFLAGS) src/SlidingWindow.cpp -o src/SlidingWindow.o

# Make the unit tests for LODESTAR.
#	Note including the utils unit tests.
.PHONY: bin/unit_tests
bin/unit_tests: bin/unit_tests/unit_test_procrustes bin/unit_tests/unit_test_vcf_parser bin/unit_tests/unit_test_sliding_window

# Create the Procrustes unit test.
bin/unit_tests/unit_test_procrustes: src/Procrustes.o utils/LinearAlgebra/MatrixOperations.o
	g++ unit_tests/unit_test_procrustes.cpp src/Procrustes.o utils/LinearAlgebra/MatrixOperations.o $(LFLAGS) bin/unit_tests/unit_test_procrustes

# Create the VCFParser unit test.
bin/unit_tests/unit_test_vcf_parser: src/VCFParser.o
	g++ unit_tests/unit_test_vcf_parser.cpp src/VCFParser.o $(LFLAGS) bin/unit_tests/unit_test_vcf_parser

# Create the SlidingWindow unit test.
bin/unit_tests/unit_test_sliding_window: src/VCFParser.o src/SlidingWindow.o 
	g++ unit_tests/unit_test_sliding_window.cpp src/VCFParser.o src/SlidingWindow.o $(LFLAGS) bin/unit_tests/unit_test_sliding_window

# Clean up object files.
.PHONY: clean
clean:
	rm bin/* bin/unit_tests/* src/*.o utils/LinearAlgebra/*.o