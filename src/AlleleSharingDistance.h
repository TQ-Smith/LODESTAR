
// File: AlleleSharingDistance.h
// Date: 
// Author: TQ Smith
// Purpose: 

#ifndef _ALLELE_SHARING_DISTANCE_H_
#define _ALLELE_SHARING_DISTANCE_H_

#include "HaplotypeEncoder.h"

static inline int IBS(Genotype s1, Genotype s2) {
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

#endif