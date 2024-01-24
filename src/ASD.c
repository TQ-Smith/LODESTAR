
#include "ASD.h"

#include "Matrix.h"

#include <stdlib.h>

#define IBS(leftA, rightA, leftB, rightB) (\
    (!!(leftA ^ rightA ^ leftB ^ rightB) ? 0 : (!(leftA ^ rightA ^ leftB ^ rightB) ? 2 : 1)))

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

    int N = asd -> numSamples;
    int index_ij;
    int ibs;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            index_ij = INDEX(i, j, N);
            
            if (encoder -> leftHaplotype[i] != encoder -> numLeaves && encoder -> leftHaplotype[j] != encoder -> numLeaves) {
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
                windowASDMatrix[i][j] = windowASDMatrix[j][i] = 1 - (asd -> curWinSim[index_ij] / (2 * asd -> curWinSimCounts[index_ij]));
                curWinSim[index_ij] = nextWinSim[index_ij];
                curWinSimCounts[index_ij] = nextWinSimCounts[index_ij];
                nextWinSim[index_ij] = 0;
                nextWinSimCounts[index_ij] = 0;
            }

            if (!isSameChromosome) {
                curWinSim[index_ij] = 0;
                curWinSimCounts[index_ij] = 0;
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
    destroy_matrix(asd -> windowASDMatrix);
    free(asd);

}