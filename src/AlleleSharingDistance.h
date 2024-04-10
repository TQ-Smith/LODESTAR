
// File: AlleleSharingDistance.h
// Date: 
// Author: TQ Smith
// Purpose: 

#ifndef _ALLELE_SHARING_DISTANCE_H_
#define _ALLELE_SHARING_DISTANCE_H_

#include "HaplotypeEncoder.h"

typedef struct {
    unsigned int ibs0;
    unsigned int ibs1;
    unsigned int ibs2;
} IBS;

static inline void reset_ibs(IBS* ibs) {
    ibs -> ibs0 = 0;
    ibs -> ibs1 = 0;
    ibs -> ibs2 = 0;
}

static inline int num_shared_alleles(Genotype s1, Genotype s2) {
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

static inline double ibs_to_asd(IBS ibs) {
    return 1.0 - (ibs.ibs1 + (2.0 * ibs.ibs2)) / (2.0 * (ibs.ibs0 + ibs.ibs1 + ibs.ibs2));
}

void process_window_multi_thread(Genotype** winGeno, IBS* winAlleleCounts, IBS* globalAlleleCounts, double* asd, int numHapsInWin, bool isLastWinOnChrom, int numSamples, int STEP_SIZE);

void process_window_single_thread(Genotype** winGeno, IBS* winAlleleCounts, IBS* stepAlleleCounts, IBS* globalAlleleCounts, double* asd, int numHapsInWin, bool isFirstWinOnChrom, bool isLastWinOnChrom, int numSamples, int STEP_SIZE);

void pairwise_ibs(Genotype* genotypes, IBS* alleleCounts, int numSamples);

#endif