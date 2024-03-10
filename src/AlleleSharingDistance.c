
// File: AlleleSharingDistance.h
// Date: 5 March 2024
// Author: TQ Smith
// Purpose: Defines the processing of a haplotype and ASD calculations.

#include "AlleleSharingDistance.h"

void process_haplotype_single_thread(double* leftHaps, double* rightHaps, double** winIBS, double** offsetIBS, double** winASD, double** globalIBS, int numSamples, int STEP_SIZE, int numHapsInWin, bool isSameChrom, int curHap) {
    int ibs;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            
        }
    }
}

void process_haplotype_multi_thread(double* leftHaps, double* rightHaps, double** winIBS, double** globalIBS, int numSamples, int STEP_SIZE, int numHapsInWin, bool isSameChrom, int curHap) {
    int ibs;
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            if (curHap == 0) {
                winIBS[i][j] = 0;
                winIBS[j][i] = 0;
            }
            if (!IS_EQUAL(leftHaps[i], MISSING) && !IS_EQUAL(leftHaps[j], MISSING)) {
                ibs = IBS(leftHaps[i], rightHaps[i], leftHaps[j], rightHaps[j]);
                winIBS[i][j] += ibs;
                winIBS[j][i]++;
                if (!isSameChrom || curHap < STEP_SIZE) {
                    globalIBS[i][j] += ibs;
                    globalIBS[j][i]++;
                }
            }
            if (curHap == numHapsInWin - 1) {
                winIBS[i][j] = 1.0 - (winIBS[i][j] / (2.0 * winIBS[j][i]));
                winIBS[j][i] = winIBS[i][j];
                winIBS[i][i] = 0;
            }
        }
    }
}