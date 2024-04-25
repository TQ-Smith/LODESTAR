
// File: AlleleSharingDistance.h
// Date: 
// Author: TQ Smith
// Purpose: 

#include "AlleleSharingDistance.h"

#include "Matrix.h"

#define increment_ibs_value(ibs, numShared) ((*(((unsigned int *) &(ibs)) + numShared))++)

void process_haplotype_multi_thread(Genotype* geno, IBS* winAlleleCounts, IBS* globalAlleleCounts, double* asd, int curHapInWin, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (curHapInWin == 0)
                reset_ibs(&winAlleleCounts[PACKED_INDEX(i, j)]);
            if (geno[i].left != MISSING && geno[j].right != MISSING) {
                numSharedAlleles = num_shared_alleles(geno[i], geno[j]);
                // printf("IBS between %d:%ld/%ld and %d:%ld/%ld is %d\n", i, geno[i].left, geno[i].right, j, geno[j].left, geno[j].right, numSharedAlleles);
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

void process_haplotype_single_thread(Genotype* geno, Genotype* stepGeno, IBS* winAlleleCounts, IBS* stepAlleleCounts, IBS* globalAlleleCounts, double* asd, int curHapInWin, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE, int WINDOW_SIZE) {
    /*
    printf("Processing Window: ");
    for (int i = 0; i < numSamples; i++)
        printf("%ld/%ld\t", geno[i].left, geno[i].right);
    printf("\n");
    printf("Processing Offset: ");
    for (int i = 0; i < numSamples; i++)
        printf("%ld/%ld\t", stepGeno[i].left, stepGeno[i].right);
    printf("\n\n");
    */
    int numSharedAlleles;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (curHapInWin == 0) {
                reset_ibs(&winAlleleCounts[PACKED_INDEX(i, j)]);
                reset_ibs(&stepAlleleCounts[PACKED_INDEX(i, j)]);
            }
            if (geno[i].left != MISSING && geno[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(geno[i], geno[j]);
                increment_ibs_value(winAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
            if (!isLastWinOnChrom && curHapInWin >= (WINDOW_SIZE - STEP_SIZE) && stepGeno[i].left != MISSING && stepGeno[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(stepGeno[i], stepGeno[j]);
                increment_ibs_value(stepAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
            if (curHapInWin == numHapsInWin - 1) {
                asd[PACKED_INDEX(i, j)] = ibs_to_asd(winAlleleCounts[PACKED_INDEX(i, j)]);
                add_ibs(&globalAlleleCounts[PACKED_INDEX(i, j)], &stepAlleleCounts[PACKED_INDEX(i, j)]);
                if (isLastWinOnChrom) {
                    add_ibs(&globalAlleleCounts[PACKED_INDEX(i, j)], &winAlleleCounts[PACKED_INDEX(i, j)]);
                } else {
                    subtract_ibs(&winAlleleCounts[PACKED_INDEX(i, j)], &stepAlleleCounts[PACKED_INDEX(i, j)]);
                    reset_ibs(&stepAlleleCounts[PACKED_INDEX(i, j)]);
                }
            }
        }
    }
}

void process_window_single_thread(Genotype** winGeno, int winStartIndex, IBS* winAlleleCounts, IBS* stepAlleleCounts, IBS* globalAlleleCounts, double* asd, int numHapsInWin, bool isFirstWinOnChrom, bool isLastWinOnChrom, int numSamples, int STEP_SIZE, int WINDOW_SIZE) {
    if (isFirstWinOnChrom) {
        int stepIndex = winStartIndex;
        for (int i = 0; i < numHapsInWin; i++) {
            process_haplotype_single_thread(winGeno[(winStartIndex + i) % WINDOW_SIZE], winGeno[stepIndex], winAlleleCounts, stepAlleleCounts, globalAlleleCounts, asd, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE, WINDOW_SIZE);
            if (i >= WINDOW_SIZE - STEP_SIZE)
                stepIndex = (stepIndex + 1) % WINDOW_SIZE;
        }
    } else {
        for (int i = WINDOW_SIZE - STEP_SIZE, j = 0; i < numHapsInWin; i++, j++)
            process_haplotype_single_thread(winGeno[(winStartIndex + i) % WINDOW_SIZE], winGeno[(winStartIndex + j) % WINDOW_SIZE], winAlleleCounts, stepAlleleCounts, globalAlleleCounts, asd, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE, WINDOW_SIZE);
    }
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