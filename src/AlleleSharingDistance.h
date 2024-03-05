
// File: AlleleSharingDistance.h
// Date: 5 March 2024
// Author: TQ Smith
// Purpose: Defines the processing of a haplotype and ASD calculations.

#ifndef _ALLELE_SHARING_DISTANCE_
#define _ALLELE_SHARING_DISTANCE_

#include <stdbool.h>

void process_haplotype_single_thread(double* leftHaps, double* rightHaps, double** winCalcs, double** overlapCalcs, double** globalCalcs, int numSamples, int numHapsInWin, bool isSameChrom);

#endif