
// File: SlidingWindow.c
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Slides a window along the contents of a VCF file.

#include "SlidingWindow.h"

#include "Matrix.h"

#include "../klib/klist.h"
#define do_not_free_window_ptr(w)
KLIST_INIT(WindowPtr, Window*, do_not_free_window_ptr)

#include <pthread.h>

pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t listLock = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    VCFGenotypeParser* parser;
    HaplotypeEncoder* encoder;
    klist_t(WindowPtr)* windowList;
    int HAP_SIZE;
    int STEP_SIZE;
    int WINDOW_SIZE;
    int curNumWin;
    int curNumWinOnChrom;
    kstring_t* curChrom;
} SharedThreadResources;

void* process_window(void* arg) {
    printf("Here from thread!\n");
    return NULL;
}

Window** window_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows) {
    
    SharedThreadResources* resources = (SharedThreadResources*) malloc(sizeof(SharedThreadResources));
    resources -> parser = parser;
    resources -> encoder = encoder;
    resources -> windowList = kl_init(WindowPtr);
    resources -> HAP_SIZE = HAP_SIZE;
    resources -> STEP_SIZE = STEP_SIZE;
    resources -> WINDOW_SIZE = WINDOW_SIZE;
    resources -> curNumWin = 1;
    resources -> curNumWinOnChrom = 1;
    resources -> curChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(ks_str(parser -> nextChromosome), resources -> curChrom);
    
    pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS - 1; i++)
        pthread_create(&threads[i], NULL, process_window, (void*) resources);
    
    printf("Here from main!\n");

    for (int i = 0; i < NUM_THREADS - 1; i++)
        pthread_join(threads[i], NULL);

    free(threads);
    kl_destroy(WindowPtr, resources -> windowList);
    free(resources -> curChrom -> s); free(resources -> curChrom);
    free(resources);
    
    return NULL;
}
