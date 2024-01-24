
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

    return asd;

}

void process_haplotype(HaplotypeEncoder* encoder, ASD* asd, int OFFSET_SIZE, bool isSameChromosome) {

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
    free(asd);

}