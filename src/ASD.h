
#ifndef _ASD_
#define _ASD_

#include "HaplotypeEncoder.h"

#include "Window.h"

#include <stdbool.h>

#define TRIANGULAR_NUMBER(N) ((N * N + N) / 2)

#define INDEX(i, j, N) ((i <= j) ? (i * N + (j - 1)) : (j * N + (i - 1)))

typedef struct {

    int numSamples;

    double* curWinSim;
    double* nextWinSim;
    double* curWinSimCounts;
    double* nextWinSimCounts;
    double* globalSim;
    double* globalSimCounts;

    double** windowASDMatrix;

} ASD;

ASD* init_asd(int numSamples);

void process_haplotype(HaplotypeEncoder* encoder, ASD* asd, int numHapsInOverlap, bool isSameChromosome, int OFFSET_SIZE, int WINDOW_SIZE);

void destroy_asd(ASD* asd);

#endif