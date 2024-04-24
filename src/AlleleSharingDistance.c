
// File: AlleleSharingDistance.h
// Date: 
// Author: TQ Smith
// Purpose: 

#include "AlleleSharingDistance.h"

#include "Matrix.h"

#define increment_ibs_value(ibs, numShared) ((*(((unsigned int *) &(ibs)) + numShared))++)

void process_haplotype_multi_thread(Genotype* genotypes, IBS* winAlleleCounts, IBS* globalAlleleCounts, double* asd, int curHapInWin, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (curHapInWin == 0)
                reset_ibs(&winAlleleCounts[PACKED_INDEX(i, j)]);
            if (genotypes[i].left != MISSING && genotypes[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(genotypes[i], genotypes[j]);
                increment_ibs_value(winAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
                if (curHapInWin < STEP_SIZE || isLastWinOnChrom)
                    increment_ibs_value(globalAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
            if (curHapInWin == numHapsInWin - 1)
                asd[PACKED_INDEX(i, j)] = ibs_to_asd(winAlleleCounts[PACKED_INDEX(i, j)]);
        }
    }
}

void process_window_multi_thread(Genotype** winGeno, IBS* winAlleleCounts, IBS* globalAlleleCounts, double* asd, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    for (int i = 0; i < numHapsInWin; i++)
        process_haplotype_multi_thread(winGeno[i], winAlleleCounts, globalAlleleCounts, asd, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE);
}

void process_haplotype_single_thread() {
    
}

void process_window_single_thread() {
    
}

void pairwise_ibs(Genotype* geno, IBS* alleleCounts, int numSamples) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) 
        for (int j = i + 1; j < numSamples; j++) 
            if (geno[i].left != MISSING && geno[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(geno[i], geno[j]);
                increment_ibs_value(alleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
}