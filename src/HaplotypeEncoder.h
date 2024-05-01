
// File: HaplotypeEncoder.h
// Date: 
// Author: TQ Smith
// Purpose: 

#ifndef _HAPLOTYPE_ENCODER_H_
#define _HAPLOTYPE_ENCODER_H_

#include <stdlib.h>

#include <stdbool.h>

#include "VCFLocusParser.h"

// Initialize klib hash table.
#include "../lib/khash.h"
KHASH_MAP_INIT_INT(haplotype, unsigned long)

#define MISSING 0xFFFFFFFFFFFFFFFF

typedef unsigned long Haplotype;

typedef struct {
    Haplotype left;
    Haplotype right;
} Genotype;

typedef struct {

    int numSamples;
    Locus* locus;
    Genotype* genotypes;

    kstring_t* chrom;
    unsigned int startLocus;
    unsigned int endLocus;
    int numLoci;

    khash_t(haplotype)* labelMap;

    unsigned long numLeaves;

} HaplotypeEncoder;

HaplotypeEncoder* init_haplotype_encoder(int numSamples);

bool get_next_haplotype(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE);

void destroy_haplotype_encoder(HaplotypeEncoder* encoder);

#endif