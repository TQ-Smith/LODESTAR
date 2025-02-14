
// File: HaplotypeEncoder.h
// Date: 6 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Assume a VCF file is phased and encode haplotypes of samples' genotypes.

#ifndef _HAPLOTYPE_ENCODER_H_
#define _HAPLOTYPE_ENCODER_H_

#include "VCFLocusParser.h"

// Initialize klib hash table.
//  Used to relabel haplotype encodings.
#include "../lib/khash.h"
KHASH_MAP_INIT_INT(haplotype, unsigned long)

// A haplotype is encoded as an unsigned long.
typedef unsigned long Haplotype;

// We collapse a haplotype with missing data to the same haplotype encoding.
#define MISSING 0xFFFFFFFFFFFFFFFF

// We define the genotype of a sample as the encodings of their left and right haplotype.
typedef struct {
    Haplotype left;
    Haplotype right;
} Genotype_t;

// Encoding scheme:
//  Right now, the haplotypes are encoded using a simplified version
//  of arithmetic coding. We form a hypothetical tree where a path from
//  the root of the tree to a leaf defines a haplotype. The tree grows
//  exponentially with the number of genotypes read in. Once the tree 
//  grows too large, we relabel the leaves starting from label 0 to
//  the number of unique haplotypes. There are 2N possible unique haplotypes,
//  where N is the number of samples, which is much less than 2^64.
//  All haplotypes with missing genotypes are treated as MISSING and ignored
//  in the labeling. This is fast for short haplotypes where each locus has
//  a small amount of alleles.

typedef struct {
    // Number of samples whose haplotypes we are tracking.
    int numSamples;
    // A locus auxilary array to read in genotypes from a VCFLocusParser.
    //  Make use of HaplotypeEncoder_t easier to work with.
    Locus* locus;
    Genotype_t* genotypes;

    // The chromsome of the current haplotype.
    kstring_t* chrom;
    // The start and end coordinates of the current haplotype.
    unsigned int startCoord;
    unsigned int endCoord;
    // The number of loci in the current haplotype.
    int numLoci;

    // A hash table used to relabel the haplotypes.
    khash_t(haplotype)* labelMap;

    // The number of possible haplotypes based off the genotypes 
    //  read in for the current haplotype. 
    unsigned long numLeaves;

} HaplotypeEncoder_t;

// Create a haplotype encoder.
// Accepts:
//  int numSamples -> The number of samples to track.
// Returns: HaplotypeEncoder_t*, The created encoder.
HaplotypeEncoder_t* init_haplotype_encoder(int numSamples);

// Encode the next haplotype from the parser.
// Accepts:
//  VCFLocusParser_t* parser -> The parser we are using the read in the VCF file.
//  HaplotypeEncoder_t* encoder -> The encoder we are encoding the haplotypes with.
//  int HAP_SIZE -> The maximum number of loci that can be encoded per haplotype.
// Returns: bool, True when the parser is not EOF and there exists another haplotype on
//              the same chromosome as the current haplotype. Otherwise, false. This makes
//              the windowing operation cleaner.
bool get_next_haplotype(VCFLocusParser_t* parser, HaplotypeEncoder_t* encoder, int HAP_SIZE);

// Destroy a haplotype encoder.
// Accepts:
//  HaplotypeEncoder_t* encoder -> The encoder to destroy.
// Returns: void.
void destroy_haplotype_encoder(HaplotypeEncoder_t* encoder);

#endif