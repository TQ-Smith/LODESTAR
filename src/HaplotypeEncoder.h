
// File: HaplotypeEncoder.h
// Date: 18 Janurary 2023
// Author: TQ Smith
// Purpose: Track the haplotypes of each sample using a simplified arithmetic encoding.

#ifndef _HAPLOTYPE_TREE_
#define _HAPLOTYPE_TREE_

#include <stdlib.h>

#include <stdbool.h>

#include "VCFGenotypeParser.h"

// Initialize klib hash table.
#include "../klib/khash.h"
KHASH_MAP_INIT_INT(32, int)

// Max number of possible haplotypes before
//  the algorithm will prune and relabel the tree.
#define MAX_NUM_LEAVES (1 << 25)

// A structure to represent the encoder.
typedef struct {

    // The number of samples to encode.
    int numSamples;
    // An array to hold the genotypes read in by the VCF parser.
    //  Not explicitly needed, but it helps usage.
    GENOTYPE* genotypes;
    // Arrays to hold the left and right haplotype encodings for each sample.
    unsigned int* leftHaplotype;
    unsigned int* rightHaplotype;

    // Number of loci in the haplotye.
    int numLoci;
    // Chromosome the haplotype rests on.
    kstring_t* chromosome;
    // The start locus of the haplotype.
    int startLocus;
    // The end locus of the haplotype.
    int endLocus;

    // A hash table used to relabel the encodings.
    khash_t(32)* labelMap;

    // The number of leaves in the haplotype tree.
    int numLeaves;

} HaplotypeEncoder;

// Creates a HaplotypeEncoder structure.
// Accepts:
//  int numSamples -> The number of samples to track.
// Returns:
//  HaplotypeEncoder*, The created structure.
HaplotypeEncoder* init_haplotype_encoder(int numSamples);

// Read in the next haplotype from a VCF file.
// Accepts:
//  VCFGenotypeParser* parser -> The parser for the VCF file.
//  HaplotypeEncoder* encoder -> The HaplotypeEncoder used to label unique haplotypes.
//  bool collapseMissingGenotypes -> Is set, all samples with a missing genotype are set to the same haplotype, which is numLeaves.
//  int HAP_SIZE -> The maximum number of loci that defines a haplotype.
// Returns:
//  bool, Returns true if the haplotype contains HAP_SIZE loci and the next loci is on the same chromosome and EOF was not reached.
bool get_next_haplotype(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, bool collapseMissingGenotypes, int HAP_SIZE);

// Relabels haplotype encodings. Simplifies tree.
// Accepts:
//  HaplotypeEncoder* encoder -> The encoder to simplify.
// Returns:
//  void.
void relabel_haplotypes(HaplotypeEncoder* encoder);

// Deallocated memory used by the encoder.
// Accepts:
//  HaplotypeEncoder* encoder -> The encoder to deallocate.
// Returns:
//  void.
void destroy_haplotype_encoder(HaplotypeEncoder* encoder);

#endif