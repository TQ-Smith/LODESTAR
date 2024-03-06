
// File: AlleleSharingDistance.h
// Date: 5 March 2024
// Author: TQ Smith
// Purpose: Defines the processing of a haplotype and ASD calculations.

#include "AlleleSharingDistance.h"

void process_haplotype_single_thread(HaplotypeEncoder* encoder, double** winIBS, double** overlapIBS, double** globalIBS, double** asdCalcs, int numHapsInWin, bool isSameChrom, int STEP_SIZE, int WINDOW_SIZE) {
    int ibs;
    for (int i = 0; i < encoder -> numSamples; i++) {
        for (int j = 0; j < encoder -> numSamples; j++) {
            if (encoder -> leftHaps[i] != MISSING && encoder -> leftHaps[j] != MISSING) {
                ibs = IBS(encoder -> leftHaps[i], encoder -> rightHaps[i], encoder -> leftHaps[j], encoder -> rightHaps[j]);
                winIBS[i][j] += ibs;
                winIBS[j][i]++;
                globalIBS[i][j] += ibs;
                globalIBS[j][i]++;
                if (numHapsInWin >= STEP_SIZE) {
                    overlapIBS[i][j] += ibs;
                    overlapIBS[j][i]++;
                }
            }
            if (!isSameChrom || numHapsInWin == WINDOW_SIZE - 1) {
                asdCalcs[i][j] = asdCalcs[j][i] = 1 - (winIBS[i][j] / (2 * winIBS[j][i]));
                asdCalcs[i][i] = 0;
                winIBS[i][j] = overlapIBS[i][j];
                winIBS[j][i] = overlapIBS[j][i];
                overlapIBS[i][j] = overlapIBS[j][i] = 0;
                if (!isSameChrom)
                    winIBS[i][j] = winIBS[j][i] = 0;
            }
        }
    }
}