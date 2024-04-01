
// File: AlleleSharingDistance.h
// Date: 
// Author: TQ Smith
// Purpose: 

#include "AlleleSharingDistance.h"

void pairwise_ibs(IBS** ibs, Genotype* genotypes, int numSamples) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (genotypes[i].left != MISSING && genotypes[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(genotypes[i], genotypes[j]);
                (*(((unsigned int *) &(ibs[i][j])) + numSharedAlleles))++;
            }
        }
    }
}