
// File: HaplotypeEncoder.h
// Date: 
// Author: TQ Smith
// Purpose: 

#ifndef _HAPLOTYPE_TREE_
#define _HAPLOTYPE_TREE_

#include <stdlib.h>

#include <stdbool.h>

#include "VCFGenotypeParser.h"

#define MAX_NUMBER_OF_LOCI 100

#define MISSING (-1)

typedef struct {

    int numSamples;

    GENOTYPE* genos;

    double* leftHaps;
    double* rightHaps;

    int numLoci;

    kstring_t* chrom;

    int startLocus;

    int endLocus;

} HaplotypeEncoder;


HaplotypeEncoder* init_haplotype_encoder(int numSamples);


bool get_next_haplotype(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, bool collapseMissingGenotypes, int HAP_SIZE);


void destroy_haplotype_encoder(HaplotypeEncoder* encoder);

#endif