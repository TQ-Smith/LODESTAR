
// File: SlidingWindow.c
// Date: 
// Author: TQ Smith
// Purpose: 

#include "SlidingWindow.h"

#include "../klib/klist.h"
#define do_not_free_window_ptr(w)
KLIST_INIT(WindowPtr, Window*, do_not_free_window_ptr)

#include <pthread.h>

#include "AlleleSharingDistance.h"

#include "Matrix.h"
MATRIX_INIT(doubles, double)
MATRIX_INIT(ibs, IBS)

#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t windowListLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

Window* get_next_window() {
    return NULL;
}

void* sliding_window_multi_thread(void* arg) {
    return NULL;
}

void sliding_window_single_thread() {

}

Window** sliding_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows) {
    return NULL;
}

typedef struct {
    VCFLocusParser* parser;
    HaplotypeEncoder* encoder;
    IBS** ibs;
    int numLoci;
    int HAP_SIZE;
} GlobalWindowRecord;

void* global_window_multi_thread(void* arg) {
    GlobalWindowRecord* record = (GlobalWindowRecord*) arg;
    int numSamples = record -> encoder -> numSamples;

    IBS** ibs = create_ibs_matrix(numSamples, numSamples);
    Genotype* genotypes = (Genotype*) calloc(numSamples, sizeof(Genotype));
    Genotype* temp;

    pthread_mutex_lock(&genomeLock);
    while (!record -> parser -> isEOF) {
        get_next_haplotype(record -> parser, record -> encoder, true, record -> HAP_SIZE);
        record -> numLoci += record -> encoder -> numLoci;
        SWAP(genotypes, record -> encoder -> genotypes, temp);
        pthread_mutex_unlock(&genomeLock);

        pairwise_ibs(ibs, genotypes, numSamples);
    }
    pthread_mutex_unlock(&genomeLock);

    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            record -> ibs[i][j].ibs0 += ibs[i][j].ibs0;
            record -> ibs[i][j].ibs1 += ibs[i][j].ibs1;
            record -> ibs[i][j].ibs2 += ibs[i][j].ibs2;
        }
    }
    pthread_mutex_unlock(&globalLock);

    destroy_ibs_matrix(ibs, numSamples);
    free(genotypes);
    return NULL;
}

Window* global_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int NUM_THREADS) {

    IBS** ibs = create_ibs_matrix(encoder -> numSamples, encoder -> numSamples);
    double** asd = create_doubles_matrix(encoder -> numSamples, encoder -> numSamples);

    Window* window = init_window();
    window -> winNum = 0;
    window -> winNumOnChrom = 0;
    kputs("Global", window -> chromosome);
    window -> startLocus = 0;
    window -> endLocus = 0;

    if (NUM_THREADS == 1) {
        while (!parser -> isEOF) {
            get_next_haplotype(parser, encoder, true, HAP_SIZE);
            window -> numLoci += encoder -> numLoci;
            pairwise_ibs(ibs, encoder -> genotypes, encoder -> numSamples);
        }
    } else {
        GlobalWindowRecord* record = (GlobalWindowRecord*) calloc(1, sizeof(GlobalWindowRecord));
        record -> parser = parser;
        record -> encoder = encoder;
        record -> ibs = ibs;
        record -> numLoci = 0;
        record -> HAP_SIZE = HAP_SIZE;

        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, global_window_multi_thread, (void*) record);
    
        global_window_multi_thread((void*) record);

        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);

        window -> numLoci = record -> numLoci;

        free(threads);
        free(record);
    }

    printf("Global ASD:\n");
    for (int i = 0; i < encoder -> numSamples; i++) {
        for (int j = i + 1; j < encoder -> numSamples; j++) {
            printf("%5lf\t", ibs_to_asd(ibs[i][j]));
        }
        printf("\n");
    }

    destroy_ibs_matrix(ibs, encoder -> numSamples);
    destroy_doubles_matrix(asd, encoder -> numSamples);

    return window;

}