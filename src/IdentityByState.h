
// File: IdentityByState.h
// Date: 8 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines allele identity-by-state counts between samples.

#ifndef _IDENTITY_BY_STATE_H_
#define _IDENTITY_BY_STATE_H_

#include <math.h>

// Our structure that tracks the counts of loci
//  with 0, 1, and 2 alleles IBS between samples.
typedef struct {
    unsigned int ibs0;
    unsigned int ibs1;
    unsigned int ibs2;
} IBS_t;

// Sets IBS counts to 0.
// Accepts:
//  IBS_t* ibs -> Sets counts to 0.
// Returns: void.
static inline void reset_ibs(IBS_t* ibs) {
    ibs -> ibs0 = 0;
    ibs -> ibs1 = 0;
    ibs -> ibs2 = 0;
}

// Adds counts from right to left.
// Accepts:
//  IBS_t* left -> The accumulating counts.
//  IBS_t* right -> The increment counts.
// Returns: void.
static inline void add_ibs(IBS_t* left, IBS_t* right) {
    left -> ibs0 += right -> ibs0;
    left -> ibs1 += right -> ibs1;
    left -> ibs2 += right -> ibs2;
}

// Subtracts counts from left by right.
// Accepts:
//  IBS_t* left -> The accumulating counts.
//  IBS_t* right -> The decrement counts.
// Returns: void.
static inline void subtract_ibs(IBS_t* left, IBS_t* right) {
    left -> ibs0 -= right -> ibs0;
    left -> ibs1 -= right -> ibs1;
    left -> ibs2 -= right -> ibs2;
}

// Converts IBS counts to asd.
// Accepts:
//  IBS_t ibs -> The IBS counts.
// Returns: double, The ln of asd represented by the IBS counts.
//              -1, if ASD is undefined.
static inline double ibs_to_asd(IBS_t ibs) {
    if (ibs.ibs0 == 0 && ibs.ibs1 == 0 && ibs.ibs2 == 0)
        return -1;
    double asd = 1.0 - (ibs.ibs1 + (2.0 * ibs.ibs2)) / (2.0 * (ibs.ibs0 + ibs.ibs1 + ibs.ibs2));
    return asd;
    // if (asd == 0)
    //    return 1;
    //else 
    //    return -log(asd);
    //return (ibs.ibs1 + (2.0 * ibs.ibs2));
}

#endif