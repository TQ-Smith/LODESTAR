#
# File: Makefile
# Date: 22 July 2023
# Author: TQ Smith
# Principal Investigator: Dr. Zachary Szpiech
# Purpose: Make and run all unit tests for the 'utils' library.
#
# Note: This could be made more eloquent by removing the ../LinearAlgebra prefix with
#			some variable creation.
#

CFLAGS = -std=c++11 -c -Wall -g
LFLAGS = -g -o

# Compile all tests, run tests, and then clean up object files and temporary .dSYM files.
.PHONY: all
all: unit_tests run_tests

# Compile all of the unit tests.
.PHONY: unit_tests
unit_tests: unit_test_command_line_argument_parser unit_test_matrix_operations unit_test_multidimensional_scaling

# Make the command line argument parser unit tests.
unit_test_command_line_argument_parser:
	g++ $(LFLAGS) bin/unit_test_command_line_argument_parser unit_test_command_line_argument_parser.cpp

# Make the matrix operations unit tests.
unit_test_matrix_operations: ../LinearAlgebra/MatrixOperations.o unit_test_matrix_operations.o
	g++ $(LFLAGS) bin/unit_test_matrix_operations unit_test_matrix_operations.o ../LinearAlgebra/MatrixOperations.o

# Compile unit tests for matrix operations.
unit_test_matrix_operations.o:
	g++ $(CFLAGS) unit_test_matrix_operations.cpp -o unit_test_matrix_operations.o

# Make the multidimensional scaling unit tests.
unit_test_multidimensional_scaling: ../LinearAlgebra/MatrixOperations.o ../LinearAlgebra/MultidimensionalScaling.o ../LinearAlgebra/NumericalRecipesInC.o unit_test_multidimensional_scaling.o
	g++ $(LFLAGS) bin/unit_test_multidimensional_scaling unit_test_multidimensional_scaling.o ../LinearAlgebra/*.o 

# Compile unit tests for multidimensional scaling.
unit_test_multidimensional_scaling.o:
	g++ $(CFLAGS) unit_test_multidimensional_scaling.cpp -o unit_test_multidimensional_scaling.o

# Compile the matrix operations source.
../LinearAlgebra/MatrixOperations.o: 
	g++ $(CFLAGS) ../LinearAlgebra/MatrixOperations.cpp -o ../LinearAlgebra/MatrixOperations.o

# Compile the multidimensional scaling source.
../LinearAlgebra/MultidimensionalScaling.o: 
	g++ $(CFLAGS) ../LinearAlgebra/MultidimensionalScaling.cpp -o ../LinearAlgebra/MultidimensionalScaling.o

# Compile the numerical recipes in C source.
../LinearAlgebra/NumericalRecipesInC.o: 
	g++ $(CFLAGS) ../LinearAlgebra/NumericalRecipesInC.cpp -o ../LinearAlgebra/NumericalRecipesInC.o

# Run all of the unit test programs.
.PHONY: run_tests
run_tests: run_unit_test_command_line_argument_parser run_unit_test_matrix_operations run_unit_test_multidimensional_scaling

# Run the command line argument parser unit tests.
.PHONY: run_unit_test_command_line_argument_parser
run_unit_test_command_line_argument_parser:
	@echo
	@echo "*****Running Unit Tests for CommandLineArgumentParser*****"
	./bin/unit_test_command_line_argument_parser
	@echo "*****Finished Unit Tests for CommandLineArgumentParser*****"

# Run the matric operations unit tests.
.PHONY: run_unit_test_matrix_operations
run_unit_test_matrix_operations:
	@echo
	@echo "*****Running Unit Tests for MatrixOperations*****"
	./bin/unit_test_matrix_operations
	@echo "*****Finished Unit Tests for MatrixOperations*****"

# Run the multidimensional scaling unit tests.
.PHONY: run_unit_test_multidimensional_scaling
run_unit_test_multidimensional_scaling:
	@echo
	@echo "*****Running Unit Tests for MultidimensionalScaling*****"
	./bin/unit_test_multidimensional_scaling
	@echo "*****Finished Unit Tests for MultidimensionalScaling*****"

# Clean everything up. Grrr...Clang
.PHONY: clean
clean:
	rm -rf ../LinearAlgebra/*.o *.dSYM bin/* *.o
