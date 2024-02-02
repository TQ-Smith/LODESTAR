
#include "ASD.h"

#include "Matrix.h"

#include <stdlib.h>

ASD* init_asd(int numSamples) {

    ASD* asd = (ASD*) calloc(1, sizeof(ASD));

    asd -> numSamples = numSamples;

    asd -> curWinSim = (double*) calloc(TRIANGULAR_NUMBER(numSamples), sizeof(double));
    asd -> nextWinSim = (double*) calloc(TRIANGULAR_NUMBER(numSamples), sizeof(double));
    asd -> curWinSimCounts = (double*) calloc(TRIANGULAR_NUMBER(numSamples), sizeof(double));
    asd -> nextWinSimCounts = (double*) calloc(TRIANGULAR_NUMBER(numSamples), sizeof(double));
    asd -> globalSim = (double*) calloc(TRIANGULAR_NUMBER(numSamples), sizeof(double));
    asd -> globalSimCounts = (double*) calloc(TRIANGULAR_NUMBER(numSamples), sizeof(double));

    asd -> windowASDMatrix = create_matrix(numSamples, numSamples);

    return asd;

}

void process_haplotype(HaplotypeEncoder* encoder, ASD* asd, int numHapsInOverlap, bool isSameChromosome, int OFFSET_SIZE, int WINDOW_SIZE) {

    int index_ij, ibs;
    for (int i = 0; i < asd -> numSamples; i++) {
        for (int j = i + 1; j < asd -> numSamples; j++) {
            index_ij = INDEX(i, j, asd -> numSamples);
            if (encoder -> leftHaplotype[i] != encoder -> numLeaves - 1 && encoder -> leftHaplotype[j] != encoder -> numLeaves - 1) {
                ibs = IBS(encoder -> leftHaplotype[i], encoder -> rightHaplotype[i], encoder -> leftHaplotype[j], encoder -> rightHaplotype[j]);
                asd -> curWinSim[index_ij] += ibs;
                asd -> curWinSimCounts[index_ij]++;
                asd -> globalSim[index_ij] += ibs;
                asd -> globalSimCounts[index_ij]++;
                if (numHapsInOverlap >= OFFSET_SIZE) {
                    asd -> nextWinSim[index_ij] += ibs;
                    asd -> nextWinSimCounts[index_ij]++;
                }
            }
            if (!isSameChromosome || numHapsInOverlap == WINDOW_SIZE - 1) {
                asd -> windowASDMatrix[i][j] = asd -> windowASDMatrix[j][i] = 1 - (asd -> curWinSim[index_ij] / (2 * asd -> curWinSimCounts[index_ij]));
                asd -> curWinSim[index_ij] = asd -> nextWinSim[index_ij];
                asd -> curWinSimCounts[index_ij] = asd -> nextWinSimCounts[index_ij];
                asd -> nextWinSim[index_ij] = 0;
                asd -> nextWinSimCounts[index_ij] = 0;
            }
            if (!isSameChromosome) {
                asd -> curWinSim[index_ij] = 0;
                asd -> curWinSimCounts[index_ij] = 0;
            }
        }
    }

}

void destroy_asd(ASD* asd) {

    if (asd == NULL)
        return;
    
    free(asd -> curWinSim);
    free(asd -> nextWinSim);
    free(asd -> curWinSimCounts);
    free(asd -> nextWinSimCounts);
    free(asd -> globalSim);
    free(asd -> globalSimCounts);
    destroy_matrix(asd -> windowASDMatrix, asd -> numSamples);
    free(asd);

}