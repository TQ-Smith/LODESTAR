
// File: AlleleSharingDistance.c
// Date: 8 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Calculate pairwise ASD values for genome-wide and local windows.

#include "AlleleSharingDistance.h"
// We could include Matrix.h, but we are only using this definition.
#define PACKED_INDEX(i, j) (i + j * (j + 1) / 2)

// Process a single haplotype when using multiple threads.
// Accepts:
//  Genotype_t* geno -> The genotypes at a haplotype.
//  IBS_t* winAlleleCounts -> Current window IBS counts.
//  IBS_t* globalAlleleCounts -> Genome-wide IBS counts.
//  double* asd -> Packed matrix for ASD calculations.
//  Window_t* window -> The current window. Passed to save IBS counts if flag set.
//  int curHapInWin -> The haplotype number in the window.
//  int numHapsInWin -> The number of haplotypes in the current window.
//  bool isLastWinOnChrom -> Set if the current window is the last one on the chromosome.
//  int numSamples -> The number of samples.
//  int STEP_SIZE -> The number of haplotypes in the step.
// Returns: void. 
void process_haplotype_multi_thread(Genotype_t* geno, IBS_t* winAlleleCounts, IBS_t* globalAlleleCounts, double* asd, Window_t* window, int curHapInWin, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    int numSharedAlleles;
    // For each pair ...
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            // If this is the first haplotype in the window, we set counts to 0.
            if (curHapInWin == 0) {
                reset_ibs(&winAlleleCounts[PACKED_INDEX(i, j)]);
            }
            // If neither of the genotypes are missing, we calculate IBS.
            if (geno[i].left != MISSING && geno[j].right != MISSING) {
                numSharedAlleles = num_shared_alleles(geno[i], geno[j]);
                increment_ibs_value(&winAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
                // If the haplotype is not in the overlap or it is the last window on the chromosome,
                //  we add IBS count to the global matrix.
                if (curHapInWin < STEP_SIZE || isLastWinOnChrom)
                    increment_ibs_value(&globalAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
            // If the haplotype is the last one in the window, convert counts to ASD.
            if (curHapInWin == numHapsInWin - 1) {
                if (window -> saveIBS)
                    window -> ibs[PACKED_INDEX(i, j)] = winAlleleCounts[PACKED_INDEX(i, j)];
                asd[PACKED_INDEX(i, j)] = ibs_to_asd(winAlleleCounts[PACKED_INDEX(i, j)]);
            }
        }
    }
}

void process_window_multi_thread(Genotype_t** winGeno, IBS_t* winAlleleCounts, IBS_t* globalAlleleCounts, double* asd, Window_t* window, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE) {
    // For each of the haplotypes in the window, calculate IBS counts.
    for (int i = 0; i < numHapsInWin; i++)
        process_haplotype_multi_thread(winGeno[i], winAlleleCounts, globalAlleleCounts, asd, window, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE);
}

// Process a single haplotype while using one thread.
// Accepts:
//  Genotype_t* geno -> The genotypes at the current haplotype.
//  Genotype_t* stepGeno -> The genotypes in the step region of the current window.
//                              Not referenced if first window on chromosome and 
//                              curHapInWin < WINDOW_SIZE - STEP_SIZE.
//  int winStartIndex -> The index of the first genotype in winGeno.
//  IBS_t* winAlleleCounts -> Current window IBS counts.
//  IBS_t* stepAlleleCounts -> The IBS counts in the step region.
//  IBS_t* globalAlleleCounts -> Genome-wide IBS counts.
//  double* asd -> Packed matrix for ASD calculations.
//  Window_t* window -> The current window. Passed to save IBS counts if flag set.
//  int curHapInWin -> The haplotype number in the window.
//  int numHapsInWin -> The number of haplotypes in the current window.
//  bool isFirstWinOnChrom -> Set if the current window is the first one on the chromosome.
//  bool isLastWinOnChrom -> Set if the current window is the last one on the chromosome.
//  int numSamples -> The number of samples.
//  int STEP_SIZE -> The number of haplotypes in the step.
//  int WINDOW_SIZE -> The number of haplotypes in the window.
// Return: void.
void process_haplotype_single_thread(Genotype_t* geno, Genotype_t* stepGeno, IBS_t* winAlleleCounts, IBS_t* stepAlleleCounts, IBS_t* globalAlleleCounts, double* asd, Window_t* window, int curHapInWin, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE, int WINDOW_SIZE) {
    int numSharedAlleles;
    // For each pair ...
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            // If the haplotype is the first in the window, we set counts to 0.
            if (curHapInWin == 0) {
                reset_ibs(&winAlleleCounts[PACKED_INDEX(i, j)]);
                reset_ibs(&stepAlleleCounts[PACKED_INDEX(i, j)]);
            }
            // If the haplotypes are not missing, calculate IBS and increment.
            if (geno[i].left != MISSING && geno[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(geno[i], geno[j]);
                increment_ibs_value(&winAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
            // Calculate IBS counts in step region of current window.
            if (!isLastWinOnChrom && curHapInWin >= (WINDOW_SIZE - STEP_SIZE) && stepGeno[i].left != MISSING && stepGeno[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(stepGeno[i], stepGeno[j]);
                increment_ibs_value(&stepAlleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
            // If we reached the end of the current window ...
            if (curHapInWin == numHapsInWin - 1) {
                if (window -> saveIBS)
                    window -> ibs[PACKED_INDEX(i, j)] = winAlleleCounts[PACKED_INDEX(i, j)];
                // Convert to ASD.
                asd[PACKED_INDEX(i, j)] = ibs_to_asd(winAlleleCounts[PACKED_INDEX(i, j)]);
                add_ibs(&globalAlleleCounts[PACKED_INDEX(i, j)], &stepAlleleCounts[PACKED_INDEX(i, j)]);
                if (isLastWinOnChrom) {
                    add_ibs(&globalAlleleCounts[PACKED_INDEX(i, j)], &winAlleleCounts[PACKED_INDEX(i, j)]);
                } else {
                    // Subtract step counts from current counts to get the next window.
                    subtract_ibs(&winAlleleCounts[PACKED_INDEX(i, j)], &stepAlleleCounts[PACKED_INDEX(i, j)]);
                    reset_ibs(&stepAlleleCounts[PACKED_INDEX(i, j)]);
                }
            }
        }
    }
}

// Single threaded window processing is not as straightforward as multithreaded window processing.
//  Mainly, the difference arises because, with multiple threads, we end up recalculating some 
//  haplotype IBS counts per the number of overlapping windows. We want to avoid this while using
//  a single thread. If the window is the first on the chromosome, we consume the whole window.
//  We keep track of the IBS counts for 1 <= haplotype <= STEP_SIZE. Then, for each window after
//  the first, we subtract the IBS counts in the step region, process WINDOW_SIZE - STEP_SIZE additional
//  haplotypes for the current window, and while doing so, calculate the IBS counts for 1 <= haplotype <= STEP_SIZE
//  of the current window. We are performing each IBS calculation twice, but we are keeping the number of pairwise
//  operations linear in the number of haplotypes.

void process_window_single_thread(Genotype_t** winGeno, int winStartIndex, IBS_t* winAlleleCounts, IBS_t* stepAlleleCounts, IBS_t* globalAlleleCounts, double* asd, Window_t* window, int numHapsInWin, bool isFirstWinOnChrom, bool isLastWinOnChrom, int numSamples, int STEP_SIZE, int WINDOW_SIZE) {
    // If the first window on the chromosome, we process the whole window and save IBS counts in the step region.
    if (isFirstWinOnChrom) {
        int stepIndex = winStartIndex;
        for (int i = 0; i < numHapsInWin; i++) {
            process_haplotype_single_thread(winGeno[(winStartIndex + i) % WINDOW_SIZE], winGeno[stepIndex], winAlleleCounts, stepAlleleCounts, globalAlleleCounts, asd, window, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE, WINDOW_SIZE);
            if (i >= (WINDOW_SIZE - STEP_SIZE))
                stepIndex = (stepIndex + 1) % WINDOW_SIZE;
        }
    // Otherwise, we advance the window and calculate stepAlleleCounts for the next window.
    } else {
        for (int i = WINDOW_SIZE - STEP_SIZE, j = 0; i < numHapsInWin; i++, j++)
            process_haplotype_single_thread(winGeno[(winStartIndex + i) % WINDOW_SIZE], winGeno[(winStartIndex + j) % WINDOW_SIZE], winAlleleCounts, stepAlleleCounts, globalAlleleCounts, asd, window, i, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE, WINDOW_SIZE);
    }
}

void pairwise_ibs(Genotype_t* geno, IBS_t* alleleCounts, int numSamples) {
    int numSharedAlleles;
    // For each pair of samples, if both do not have missing haplotypes,
    //  calculate IBS count and increment.
    for (int i = 0; i < numSamples; i++) 
        for (int j = i + 1; j < numSamples; j++) 
            if (geno[i].left != MISSING && geno[j].left != MISSING) {
                numSharedAlleles = num_shared_alleles(geno[i], geno[j]);
                increment_ibs_value(&alleleCounts[PACKED_INDEX(i, j)], numSharedAlleles);
            }
}