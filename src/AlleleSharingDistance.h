
// File: AlleleSharingDistance.h
// Date: 8 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Calculate pairwise ASD values for genome-wide and local windows.

#ifndef _ALLELE_SHARING_DISTANCE_H_
#define _ALLELE_SHARING_DISTANCE_H_

#include "HaplotypeEncoder.h"
#include "IdentityByState.h"
#include "Window.h"

// Calculate the number of shared alleles between two genotypes.
// Accepts:
//  Genotype_t s1 -> The genotype of the first sample.
//  Genotype_t s2 -> The genotype of the second sample.
// Returns: int, The number of shared alleles.
static inline int num_shared_alleles(Genotype_t s1, Genotype_t s2) {
    if (s1.left == s2.left) {
        if (s1.right == s2.right)
            return 2;
        else
            return 1;
    } else if (s1.left == s2.right) {
        if (s1.right == s2.left) 
            return 2;
        else 
            return 1;
    } else {
        if (s1.right == s2.right)
            return 1;
        else if (s1.right == s2.left)
            return 1;
        else
            return 0;
    }
}

// I do not love how these methods are written.
//  Pairwise calculations are O(n^2). My goal was to
//  do as few pairwise operations as possible. A result
//  of this is the need for many arguments. For example,
//  we must know if the window is the last one on the chromosome
//  to add non-overlapping IBS values to the global counts.

// Process a whole window while using multiple threads.
// Accepts:
//  Genotype_t** winGeno -> The haplotypes for the whole window.
//  IBS_t* winAlleleCounts -> Current window IBS counts.
//  IBS_t* globalAlleleCounts -> Genome-wide IBS counts.
//  double* asd -> Packed matrix for ASD calculations.
//  Window_t* window -> The current window. Passed to save IBS counts if flag set.
//  int numHapsInWin -> The number of haplotypes in the current window.
//  bool isLastWinOnChrom -> Set if the current window is the last one on the chromosome.
//  int numSamples -> The number of samples.
//  int STEP_SIZE -> The number of haplotypes in the step.
// Returns: void. 
void process_window_multi_thread(Genotype_t** winGeno, IBS_t* winAlleleCounts, IBS_t* globalAlleleCounts, double* asd, Window_t* window, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE);

// Process a whole window while using a single thread.
// Accepts:
//  Genotype_t** winGeno -> The haplotypes for the whole window.
//  int winStartIndex -> The index of the first genotype in winGeno.
//  IBS_t* winAlleleCounts -> Current window IBS counts.
//  IBS_t* stepAlleleCounts -> The IBS counts in the step region.
//  IBS_t* globalAlleleCounts -> Genome-wide IBS counts.
//  double* asd -> Packed matrix for ASD calculations.
//  Window_t* window -> The current window. Passed to save IBS counts if flag set.
//  int numHapsInWin -> The number of haplotypes in the current window.
//  bool isFirstWinOnChrom -> Set if the current window is the first one on the chromosome.
//  bool isLastWinOnChrom -> Set if the current window is the last one on the chromosome.
//  int numSamples -> The number of samples.
//  int STEP_SIZE -> The number of haplotypes in the step.
//  int WINDOW_SIZE -> The number of haplotypes in the window.
// Return: void.
void process_window_single_thread(Genotype_t** winGeno, int winStartIndex, IBS_t* winAlleleCounts, IBS_t* stepAlleleCounts, IBS_t* globalAlleleCounts, double* asd, Window_t* window, int numHapsInWin, bool isFirstWinOnChrom, bool isLastWinOnChrom, int numSamples, int STEP_SIZE, int WINDOW_SIZE);

// Calculate pairwise IBS counts between samples at one haplotype.
// Accepts:
//  Genotype_t* geno -> The genotypes of each sample.
//  IBS_t* alleleCounts -> Current IBS counts.
//  int numSamples -> The number of samples.
// Returns: void.
void pairwise_ibs(Genotype_t* geno, IBS_t* alleleCounts, int numSamples);

#endif