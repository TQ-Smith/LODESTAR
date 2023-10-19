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
bin/lodestar: lib utils src utils/CommandLineArgumentParser.hpp
	g++ $(LFLAGS) bin/lodestar utils/LinearAlgebra/*.o lib/*.o src/*.o -lz

# Make the prewritten libraries.
.PHONY: lib
lib: lib/gzstream.o

# Make gzstream.
lib/gzstream.o:
	g++ $(CFLAGS) lib/gzstream.c -o lib/gzstream.o

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
src: src/lodestar.o src/Procrustes.o src/VCFParser.o src/SlidingWindow.o

# Make LODESTAR main executable.
src/lodestar.o:
	g++ $(CFLAGS) src/lodestar.cpp -o src/lodestar.o

# Make Procrustes.
src/Procrustes.o:
	g++ $(CFLAGS) src/Procrustes.cpp -o src/Procrustes.o

# Make VCFParser.
src/VCFParser.o:
	g++ $(CFLAGS) src/VCFParser.cpp -o src/VCFParser.o

# Make SlidingWindow.
src/SlidingWindow.o:
	g++ $(CFLAGS) src/SlidingWindow.cpp -o src/SlidingWindow.o

# Make the unit tests for LODESTAR.
#	Note including the utils unit tests.
.PHONY: bin/unit_tests
bin/unit_tests: bin/unit_tests/unit_test_procrustes bin/unit_tests/unit_test_vcf_parser bin/unit_tests/unit_test_sliding_window

# Create the Procrustes unit test.
bin/unit_tests/unit_test_procrustes: src/Procrustes.o utils/LinearAlgebra/MatrixOperations.o unit_tests/unit_test_procrustes.o
	g++ $(LFLAGS) bin/unit_tests/unit_test_procrustes unit_tests/unit_test_procrustes.o src/Procrustes.o utils/LinearAlgebra/MatrixOperations.o

# Compile procrustes analysis unit test.
unit_tests/unit_test_procrustes.o:
	g++ $(CFLAGS) unit_tests/unit_test_procrustes.cpp -o unit_tests/unit_test_procrustes.o

# Create the VCFParser unit test.
bin/unit_tests/unit_test_vcf_parser: lib/gzstream.o src/VCFParser.o unit_tests/unit_test_vcf_parser.o
	g++ $(LFLAGS) bin/unit_tests/unit_test_vcf_parser unit_tests/unit_test_vcf_parser.o lib/gzstream.o src/VCFParser.o -lz

# Compile vcf parser unit test.
unit_tests/unit_test_vcf_parser.o:
	g++ $(CFLAGS) unit_tests/unit_test_vcf_parser.cpp -o unit_tests/unit_test_vcf_parser.o

# Create the SlidingWindow unit test.
bin/unit_tests/unit_test_sliding_window: lib/gzstream.o src/VCFParser.o src/SlidingWindow.o unit_tests/unit_test_sliding_window.o
	g++ $(LFLAGS) bin/unit_tests/unit_test_sliding_window unit_tests/unit_test_sliding_window.o lib/gzstream.o src/VCFParser.o src/SlidingWindow.o -lz

# Compile sliding window unit test.
unit_tests/unit_test_sliding_window.o:
	g++ $(CFLAGS) unit_tests/unit_test_sliding_window.cpp -o unit_tests/unit_test_sliding_window.o

# Clean up object files.
.PHONY: clean
clean:
	rm bin/lodestar unit_tests/*.o bin/unit_tests/* lib/*.o src/*.o utils/LinearAlgebra/*.o