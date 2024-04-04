
// File: HaplotypeEncoder.h
// Date: 
// Author: TQ Smith
// Purpose: 

#ifndef _HAPLOTYPE_ENCODER_H_
#define _HAPLOTYPE_ENCODER_H_

#include <stdlib.h>

#include <stdbool.h>

#include "VCFLocusParser.h"

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
    int startLocus;
    int endLocus;
    int numLoci;

    unsigned long numLeaves;

} HaplotypeEncoder;

HaplotypeEncoder* init_haplotype_encoder(int numSamples);

bool get_next_haplotype(VCFLocusParser* parser, HaplotypeEncoder* encoder, bool collapseMissingGenotypes, int HAP_SIZE);

void destroy_haplotype_encoder(HaplotypeEncoder* encoder);

#endif