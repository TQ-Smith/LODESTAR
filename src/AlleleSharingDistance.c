
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

void process_window_multi_thread(Genotype** winGeno, IBS* winAlleleCounts, IBS* globalAlleleCounts, double* asd, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    for (int i = 0; i < numHapsInWin; i++) {
        process_haplotype_multi_thread(winGeno[i], winAlleleCounts, globalAlleleCounts, asd, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE);
    }
}

void process_haplotype_single_thread(Genotype** winGeno, IBS* winAlleleCounts, IBS* stepAlleleCounts, IBS* globalAlleleCounts, double* asd, int curHap, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            
        }
    }
}

void process_window_single_thread(Genotype** winGeno, IBS* winAlleleCounts, IBS* stepAlleleCounts, IBS* globalAlleleCounts, double* asd, int numHapsInWin,  bool isFirstWinOnChrom, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    int start;
    if (isFirstWinOnChrom)
        start = 0;
    else
        start = numHapsInWin - STEP_SIZE - 1;
    for (int i = start; i < numHapsInWin; i++) {
        process_haplotype_single_thread(winGeno, winAlleleCounts, stepAlleleCounts, globalAlleleCounts, asd, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE);
    }
}

void pairwise_ibs(Genotype* genotypes, IBS* alleleCounts, int numSamples) {
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