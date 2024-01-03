
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
    unsigned int* leftHaplotype;
    unsigned int* rightHaplotype;

    khash_t(32)* labelMap;

    int numLeaves;

} HaplotypeTree;

HaplotypeTree* init_haplotype_tree(int numSamples);

void add_locus(HaplotypeTree* tree, int numAlleles, GENOTYPE* genotypes, bool collapseMissingGenotypes);

void relabel_haplotypes(HaplotypeTree* tree);

void destroy_haplotype_tree(HaplotypeTree* tree);

#endif