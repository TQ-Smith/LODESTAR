
#include "ASD.h"

#include "Matrix.h"

#include <stdlib.h>

#include <math.h>

ASD* init_asd(int numSamples) {

    ASD* asd = (ASD*) calloc(1, sizeof(ASD));

    asd -> numSamples = numSamples;

    asd -> curWinSim = (double*) calloc(TRIANGULAR_NUMBER((numSamples - 1)), sizeof(double));
    asd -> nextWinSim = (double*) calloc(TRIANGULAR_NUMBER((numSamples - 1)), sizeof(double));
    asd -> curWinSimCounts = (double*) calloc(TRIANGULAR_NUMBER((numSamples - 1)), sizeof(double));
    asd -> nextWinSimCounts = (double*) calloc(TRIANGULAR_NUMBER((numSamples - 1)), sizeof(double));
    asd -> globalSim = (double*) calloc(TRIANGULAR_NUMBER((numSamples - 1)), sizeof(double));
    asd -> globalSimCounts = (double*) calloc(TRIANGULAR_NUMBER((numSamples - 1)), sizeof(double));

    asd -> windowASDMatrix = create_matrix(numSamples, numSamples);

    return asd;

}

typedef struct {
    HaplotypeEncoder* encoder;
    ASD* asd;
    int numHapsInOverlap;
    bool isSameChromosome;
    int OFFSET_SIZE;
    int WINDOW_SIZE;
    int startIndex;
    int partitionSize;
} ASDPartition;

ASDPartition* create_asd_partition(HaplotypeEncoder* encoder, ASD* asd, int numHapsInOverlap, bool isSameChromosome, int OFFSET_SIZE, int WINDOW_SIZE, int startIndex, int partitionSize) {
    ASDPartition* partition = malloc(sizeof(ASDPartition));
    partition -> encoder = encoder;
    partition -> asd = asd;
    partition -> numHapsInOverlap = numHapsInOverlap;
    partition -> isSameChromosome = isSameChromosome;
    partition -> OFFSET_SIZE = OFFSET_SIZE;
    partition -> WINDOW_SIZE = WINDOW_SIZE;
    partition -> startIndex = startIndex;
    partition -> partitionSize = partitionSize;
    return partition;
}

void compute_asd_partition(void* arg) {

    ASDPartition* partition = (ASDPartition*) arg;
    HaplotypeEncoder* encoder = partition -> encoder;
    ASD* asd = partition -> asd;
    int numHapsInOverlap = partition -> numHapsInOverlap;
    bool isSameChromosome = partition -> isSameChromosome;
    int OFFSET_SIZE = partition -> OFFSET_SIZE;
    int WINDOW_SIZE = partition -> WINDOW_SIZE;
    int startIndex = partition -> startIndex;
    int partitionSize = partition -> partitionSize;

    int N = asd -> numSamples;
    int i = N - 2 - floor(sqrt(-8*startIndex + 4*N*(N-1)-7)/2.0 - 0.5);
    int j = startIndex + i + 1 - N*(N-1)/2 + (N-i)*((N-i)-1)/2;

    int index = startIndex, ibs;
    while (index - startIndex < partitionSize) {
        if (encoder -> leftHaplotype[i] != encoder -> numLeaves - 1 && encoder -> leftHaplotype[j] != encoder -> numLeaves - 1) {
            ibs = IBS(encoder -> leftHaplotype[i], encoder -> rightHaplotype[i], encoder -> leftHaplotype[j], encoder -> rightHaplotype[j]);
            asd -> curWinSim[index] += ibs;
            asd -> curWinSimCounts[index]++;
            asd -> globalSim[index] += ibs;
            asd -> globalSimCounts[index]++;
            if (numHapsInOverlap >= OFFSET_SIZE) {
                asd -> nextWinSim[index] += ibs;
                asd -> nextWinSimCounts[index]++;
            }
        }
        if (!isSameChromosome || numHapsInOverlap == WINDOW_SIZE - 1) {
            asd -> windowASDMatrix[i][j] = asd -> windowASDMatrix[j][i] = 1 - (asd -> curWinSim[index] / (2 * asd -> curWinSimCounts[index]));
            asd -> curWinSim[index] = asd -> nextWinSim[index];
            asd -> curWinSimCounts[index] = asd -> nextWinSimCounts[index];
            asd -> nextWinSim[index] = 0;
            asd -> nextWinSimCounts[index] = 0;
        }
        if (!isSameChromosome) {
            asd -> curWinSim[index] = 0;
            asd -> curWinSimCounts[index] = 0;
        }
        index++;
        j++;
        if (j == N) {
            i++;
            j = i + 1;
        }
    }

    free(partition);

}

void process_haplotype(HaplotypeEncoder* encoder, ASD* asd, ThreadPool_t* pool, int numHapsInOverlap, bool isSameChromosome, int OFFSET_SIZE, int WINDOW_SIZE) {

    ASDPartition* partition;

    if (pool == NULL) {
        partition = create_asd_partition(encoder, asd, numHapsInOverlap, isSameChromosome, OFFSET_SIZE, WINDOW_SIZE, 0, TRIANGULAR_NUMBER((asd -> numSamples - 1)));
        compute_asd_partition((void*) partition);
    } else {

        int numPairs = TRIANGULAR_NUMBER((asd -> numSamples - 1));
        int partitionSize = (int) (numPairs / (pool -> numThreads + 1));
        int startIndex = 0;

        for (int i = 0; i < pool -> numThreads; i++) {
            partition = create_asd_partition(encoder, asd, numHapsInOverlap, isSameChromosome, OFFSET_SIZE, WINDOW_SIZE, startIndex, partitionSize);
            thread_pool_add_work(pool, compute_asd_partition, (void*) partition);
            startIndex += partitionSize;
        }
        partition = create_asd_partition(encoder, asd, numHapsInOverlap, isSameChromosome, OFFSET_SIZE, WINDOW_SIZE, startIndex, numPairs - startIndex);
        thread_pool_add_work(pool, compute_asd_partition, (void*) partition);

        thread_pool_wait(pool);

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