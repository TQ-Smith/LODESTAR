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

// The number of elements in the upper triangle of a symmetric matrix.
#define PACKED_SIZE(N) ((N * (N + 1)) / 2)

// We store the upper triangle of a symmetric matrix in a one-dimensional array,
//  This is known as packed storage. Element Aij -> a_{i + j * (j + 1) / 2}
#define PACKED_INDEX(i, j) (i + j * (j + 1) / 2)

// Better than using an if statement. ibs0, ibs1, and ibs2 correspond the 1st, 2nd, and 3rd field
//  in the structure. We offset the address of the passed IBS structure and increment the 
//  corresponding pointer. Pointer arithemtic can cause issues on different architecture.
// #define increment_ibs_value(ibs, numShared) ((*(((unsigned int *) &(ibs)) + numShared))++)
static inline void increment_ibs_value(IBS_t* ibs, int numShared) {
    switch (numShared) {
        case 0: ibs -> ibs0++; break;
        case 1: ibs -> ibs1++; break;
        case 2: ibs -> ibs2++; break;
    }
}

// Calculate the number of shared alleles between two genotypes.
// Accepts:
//  Genotype_t s1 -> The genotype of the first sample.
//  Genotype_t s2 -> The genotype of the second sample.
// Returns: int, The number of shared alleles.
static inline int num_shared_alleles(Genotype_t s1, Genotype_t s2) {
    if (s1.left == s2.left) {
        if (s1.right == s2.right)
            return 2;
        else
            return 1;
    } else if (s1.left == s2.right) {
        if (s1.right == s2.left) 
            return 2;
        else 
            return 1;
    } else {
        if (s1.right == s2.right)
            return 1;
        else if (s1.right == s2.left)
            return 1;
        else
            return 0;
    }
}


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

// Calculate IBS in blocks along the genome.
// Accepts:
//  VCFLocusParser_t* vcfFile -> The VCF file we are reading in.
//  HaplotypeEncoder_t* encoder -> The encoder for haplotypes.
//  int numSamples -> The number of samples.
//  int blockSize -> The size of the genome block in base-pairs.
//  int haplotypeSize -> The size of the haplotypes in number of loci.
//  int NUM_THREADS -> The number of threads to use in the computation.
// Returns: BlockList_t*, the list of resulting blocks.
BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int NUM_THREADS);

// Convert IBS blocks to ASD and perform Procrustes analysis
// Accepts:
//  BlockList_t* globalList -> Blocks along the genome with precomputed pairwise IBS.
//  double** y -> Centered/normalized target matrix to perform Procrustes against. If NULL, use genome-wide matrix and do not caluclate jackknife.
//  double* y0 -> Centroid of user defined points.
//  int k -> The dimension to reduce to.
//  int dropThreshold -> If numHaps within block is less than this value, we drop it.
//  int NUM_THREADS -> The number of threads to use in the computation.
// Returns: void.
void procrustes(BlockList_t* globalList, double** y, double* y0, int k, int dropThreshold, int NUM_THREADS);

#endif