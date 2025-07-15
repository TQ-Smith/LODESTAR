// File: LODESTAR.h
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main LODESTAR operations.

#ifndef _LODESTAR_H_
#define _LODESTAR_H_

#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "BlockList.h"
#include "MatrixOperations.h"
#include <math.h>

// Converts IBS counts to asd.
// Accepts:
//  IBS_t ibs -> The IBS counts.
// Returns: double, the transformed distance.
static inline double ibs_to_asd(IBS_t ibs) {
    double L = (ibs.ibs0 + ibs.ibs1 + ibs.ibs2);
    return sqrt(4 * L - 2 * (ibs.ibs1 + 2 * ibs.ibs2));
}

// Adds counts from right to left.
// Accepts:
//  IBS_t* left -> The accumulating counts.
//  IBS_t* right -> The increment counts.
// Returns: void.
static inline void add_ibs(IBS_t* left, IBS_t* right) {
    left -> ibs0 += right -> ibs0;
    left -> ibs1 += right -> ibs1;
    left -> ibs2 += right -> ibs2;
}

BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int NUM_THREADS);

#endif