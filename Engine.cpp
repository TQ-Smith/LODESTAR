//
// File: Engine.cpp
// Started: 17 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains the main logic/flow of the analysis.
//          Also handles threading.
//

// Used to define main function operations.
#include "Engine.hpp"

// Used for allocating new matrices and deep-copies.
#include "MatrixOperations.hpp"

// Used to preform the MDS operation.
#include "MultidimensionalScaling.hpp"

// Used to preform procurstes analysis and permutation tests.
#include "Procrustes.hpp"

// Used to parse genotype data file.
#include "GenotypeFileParser.hpp"