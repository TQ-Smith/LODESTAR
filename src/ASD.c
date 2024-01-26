
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

    int ibs;
    for (int i = 0; i < asd -> numSamples; i++) {
        for (int j = i + 1; j < asd -> numSamples; j++) {

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