
// File: AlleleSharingDistance.h
// Date: 
// Author: TQ Smith
// Purpose: 

#include "AlleleSharingDistance.h"

#define increment_ibs_value(ibs, numShared) ((*(((unsigned int *) &(ibs)) + numShared))++)

void process_haplotype_multi_thread(Genotype* genotypes, IBS** alleleCounts, double** asd, int curHap, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (curHap == 0)
                reset_ibs(&alleleCounts[i][j]);
            if (genotypes[i].left != MISSING && genotypes[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(genotypes[i], genotypes[j]);
                increment_ibs_value(alleleCounts[i][j], numSharedAlleles);
                if (curHap < STEP_SIZE || isLastWinOnChrom)
                    increment_ibs_value(alleleCounts[j][i], numSharedAlleles);
            }
            if (curHap == numHapsInWin - 1) {
                asd[i][j] = asd[j][i] = ibs_to_asd(alleleCounts[i][j]);
                asd[i][i] = 0;
            }
        }
    }
}

void pairwise_ibs_single_thread() {

}


void pairwise_ibs(IBS** alleleCounts, Genotype* genotypes, int numSamples) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (genotypes[i].left != MISSING && genotypes[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(genotypes[i], genotypes[j]);
                increment_ibs_value(alleleCounts[i][j], numSharedAlleles);
            }
        }
    }
}