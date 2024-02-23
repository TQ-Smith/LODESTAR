
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

#include <unistd.h>

pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t listLock = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    VCFGenotypeParser* parser;
    HaplotypeEncoder* encoder;
    klist_t(WindowPtr)* windowList;

    int HAP_SIZE;
    int STEP_SIZE;
    int WINDOW_SIZE;

    double** overlapLeftHaps;
    double** overlapRightHaps;
    int numLociInOverlap;
    int* overlapStartLoci;

    int curWinNum;
    int curWinOnChromNum;
    kstring_t* curChrom;
} SharedThreadResources;

void* process_window(void* arg) {

    printf("Thread %u\n", (int) pthread_self());

    SharedThreadResources* resources = (SharedThreadResources*) arg;

    double** swap = NULL;
    double** overlapLeftHaps = create_matrix(resources -> WINDOW_SIZE, resources -> encoder -> numSamples);
    double** overlapRightHaps = create_matrix(resources -> WINDOW_SIZE, resources -> encoder -> numSamples);

    int numHapsInOverlap;
    bool isSameChromosome;
    Window* window;

    while (true) {

        pthread_mutex_lock(&genomeLock);
        if (resources -> parser -> isEOF) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }

        window = init_window();
        window -> windowNum = ++resources -> curWinNum;
        window -> windowNumOnChromosome = ++resources -> curWinOnChromNum;
        kputs(ks_str(resources -> curChrom), window -> chromosome);
        window -> numLoci = resources -> numLociInOverlap;

        swap = resources -> overlapLeftHaps;
        resources -> overlapLeftHaps = overlapLeftHaps;
        overlapLeftHaps = swap;
        swap = resources -> overlapRightHaps;
        resources -> overlapRightHaps = overlapRightHaps;
        overlapRightHaps = swap;
        
        numHapsInOverlap = resources -> numLociInOverlap / (resources -> HAP_SIZE);
        isSameChromosome = true;

        while(numHapsInOverlap < resources -> WINDOW_SIZE && isSameChromosome) {
            isSameChromosome = get_next_haplotype(resources -> parser, resources -> encoder, true, resources -> HAP_SIZE);

            if (numHapsInOverlap % resources -> STEP_SIZE == 0)
                resources -> overlapStartLoci[(window -> windowNumOnChromosome + numHapsInOverlap) % ((resources -> WINDOW_SIZE - resources -> STEP_SIZE) / resources -> STEP_SIZE + 1)] = resources -> encoder -> startLocus;
            window -> numLoci += resources -> encoder -> numLoci;
            numHapsInOverlap++;
        }
        window -> startLocus = resources -> overlapStartLoci[window -> windowNumOnChromosome % ((resources -> WINDOW_SIZE - resources -> STEP_SIZE) / resources -> STEP_SIZE + 1)];
        window -> endLocus = resources -> encoder -> endLocus;
        if (isSameChromosome) {
            resources -> numLociInOverlap = window -> numLoci - (resources-> HAP_SIZE * resources -> STEP_SIZE);
        } else {
            resources -> curChrom -> l = 0;
            kputs(ks_str(resources -> parser -> nextChromosome), resources -> curChrom);
            resources -> numLociInOverlap = 0;
            resources -> curWinOnChromNum = 0;
        }
        pthread_mutex_unlock(&genomeLock);

        // Process window.

        pthread_mutex_lock(&listLock);
        *kl_pushp(WindowPtr, resources -> windowList) = window;
        pthread_mutex_unlock(&listLock);

    }
    
    destroy_matrix(overlapLeftHaps, resources -> WINDOW_SIZE);
    destroy_matrix(overlapRightHaps, resources -> WINDOW_SIZE);

    return NULL;
}

void sort_windows(Window** windows, int numWindows) {
    Window* temp;
    int j;
    for (int i = 1; i < numWindows; i++) {
        temp = windows[i];
        j = i - 1;
        while (j >= 0 && windows[j] -> windowNum > temp -> windowNum) {
            windows[j + 1] = windows[j];
            j = j - 1;
        }
        windows[j + 1] = temp;
    }
}

Window** window_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows) {
    
    SharedThreadResources* resources = (SharedThreadResources*) malloc(sizeof(SharedThreadResources));
    resources -> parser = parser;
    resources -> encoder = encoder;
    resources -> windowList = kl_init(WindowPtr);
    resources -> HAP_SIZE = HAP_SIZE;
    resources -> STEP_SIZE = STEP_SIZE;
    resources -> WINDOW_SIZE = WINDOW_SIZE;
    resources -> overlapLeftHaps = create_matrix(WINDOW_SIZE, encoder -> numSamples);
    resources -> overlapRightHaps = create_matrix(WINDOW_SIZE, encoder -> numSamples);
    resources -> numLociInOverlap = 0;
    resources -> overlapStartLoci = (int*) calloc((WINDOW_SIZE - STEP_SIZE) / STEP_SIZE + 1, sizeof(int));
    resources -> curWinNum = 0;
    resources -> curWinOnChromNum = 0;
    resources -> curChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(ks_str(parser -> nextChromosome), resources -> curChrom);
    
    pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS - 1; i++)
        pthread_create(&threads[i], NULL, process_window, (void*) resources);
    
    process_window((void*) resources);

    for (int i = 0; i < NUM_THREADS - 1; i++)
        pthread_join(threads[i], NULL);

    *numWindows = resources -> curWinNum;
    Window** windows = (Window**) malloc(*numWindows * sizeof(Window*));
    int i = 0;
    for (kliter_t(WindowPtr)* it = kl_begin(resources -> windowList); it != kl_end(resources -> windowList); it = kl_next(it)) {
        windows[i] = kl_val(it);
        i++;
    }
    sort_windows(windows, *numWindows);

    free(threads);
    kl_destroy(WindowPtr, resources -> windowList);
    destroy_matrix(resources -> overlapLeftHaps, WINDOW_SIZE);
    destroy_matrix(resources -> overlapRightHaps, WINDOW_SIZE);
    free(resources -> overlapStartLoci);
    free(resources -> curChrom -> s); free(resources -> curChrom);
    free(resources);
    
    return windows;

}
