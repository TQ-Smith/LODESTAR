
// File: AlleleSharingDistance.h
// Date: 
// Author: TQ Smith
// Purpose: 

#include "AlleleSharingDistance.h"

#include "Matrix.h"

#define increment_ibs_value(ibs, numShared) ((*(((unsigned int *) &(ibs)) + numShared))++)

void process_haplotype_multi_thread(Genotype* genotypes, IBS* winAlleleCounts, IBS* globalAlleleCounts, double* asd, int curHap, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (curHap == 0)
                reset_ibs(&winAlleleCounts[PACKED_INDEX(i, j)]);
            if (genotypes[i].left != MISSING && genotypes[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(genotypes[i], genotypes[j]);
                increment_ibs_value(winAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
                if (curHap < STEP_SIZE || isLastWinOnChrom)
                    increment_ibs_value(globalAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
            if (curHap == numHapsInWin - 1) {
                asd[PACKED_INDEX(i, j)] = ibs_to_asd(winAlleleCounts[PACKED_INDEX(i, j)]);
            }
        }
    }
}

void pairwise_ibs_single_thread() {

}


void pairwise_ibs(IBS* alleleCounts, Genotype* genotypes, int numSamples) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (genotypes[i].left != MISSING && genotypes[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(genotypes[i], genotypes[j]);
                increment_ibs_value(alleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
        }
    }
}