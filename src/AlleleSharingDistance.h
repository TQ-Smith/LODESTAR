
// File: AlleleSharingDistance.h
// Date: 5 March 2024
// Author: TQ Smith
// Purpose: Defines the processing of a haplotype and ASD calculations.

#ifndef _ALLELE_SHARING_DISTANCE_
#define _ALLELE_SHARING_DISTANCE_

#include "HaplotypeEncoder.h"

#include <math.h>

#define EPS 1e-10

#define IS_EQUAL(A, B) (fabs(A - B) < EPS)

static inline int IBS(double left1, double right1, double left2, double right2) {
    if (IS_EQUAL(left1, left2)) {
        if (IS_EQUAL(right1, right2))
            return 2;
        else
            return 1;
    } else if (IS_EQUAL(left1, right2)) {
        if (IS_EQUAL(right1, left2)) 
            return 2;
        else 
            return 1;
    } else {
        if (IS_EQUAL(right1, right2))
            return 1;
        else if (IS_EQUAL(right1, left2))
            return 1;
        else
            return 0;
    }
}

void process_haplotype_single_thread(HaplotypeEncoder* encoder, double** winIBS, double** overlapIBS, double** globalIBS, double** asdCalcs, int numHapsInWin, bool isSameChrom, int STEP_SIZE, int WINDOW_SIZE);

#endif