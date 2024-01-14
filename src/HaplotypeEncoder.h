
#ifndef _HAPLOTYPE_TREE_
#define _HAPLOTYPE_TREE_

#include <stdlib.h>

#include <stdbool.h>

#include "VCFGenotypeParser.h"

#include "../klib/khash.h"
KHASH_MAP_INIT_INT(32, int)

#define MAX_NUM_LEAVES (1 << 25)

typedef struct {

    int numSamples;
    GENOTYPE* genotypes;
    unsigned int* leftHaplotype;
    unsigned int* rightHaplotype;

    int numLoci;
    kstring_t* chromosome;
    int startLocus;
    int endLocus;

    khash_t(32)* labelMap;

    int numLeaves;

} HaplotypeEncoder;

HaplotypeEncoder* init_haplotype_encoder(int numSamples);

bool get_next_haplotype(VCFGenotypeParser* parser, HaplotypeEncoder* tree, int HAP_SIZE);

void relabel_haplotypes(HaplotypeEncoder* tree);

void destroy_haplotype_encoder(HaplotypeEncoder* tree);

#endif