// File: LODESTAR.c
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main LODESTAR operations.

#include "LODESTAR.h"
#include <pthread.h>

// Index lower triangle for allele counts.
#define PACKED_INDEX(i, j) (i + j * (j + 1) / 2)

BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int NUM_THREADS) {
    
}