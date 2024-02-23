
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
    int NUM_THREADS;

    double** leftHaps;
    double** rightHaps;
    double** overlapLeftHaps;
    double** overlapRightHaps;

    double** windowCalculations;
    double** overlapCalculations;
    double** globalCalculations;

    int numLociInOverlap;
    int numLociInGlobal;
    int* overlapStartLoci;

    int curWinNum;
    int curWinOnChromNum;
    kstring_t* curChrom;
} WindowRecord;

Window* get_next_window(WindowRecord* record) {

    Window* window = init_window();
    window -> windowNum = ++record -> curWinNum;
    window -> windowNumOnChromosome = ++record -> curWinOnChromNum;
    kputs(ks_str(record -> curChrom), window -> chromosome);
    window -> numLoci = record -> numLociInOverlap;
    
    int numHapsInOverlap = record -> numLociInOverlap / (record -> HAP_SIZE);
    bool isSameChromosome = true;

    while(numHapsInOverlap < record -> WINDOW_SIZE && isSameChromosome) {
        isSameChromosome = get_next_haplotype(record -> parser, record -> encoder, true, record -> HAP_SIZE);

        if (record -> NUM_THREADS == 1) {
            // TODO.
        } else {
            // TODO.
        }

        if (numHapsInOverlap % record -> STEP_SIZE == 0)
            record -> overlapStartLoci[(window -> windowNumOnChromosome + numHapsInOverlap) % ((record -> WINDOW_SIZE - record -> STEP_SIZE) / record -> STEP_SIZE + 1)] = record -> encoder -> startLocus;
        window -> numLoci += record -> encoder -> numLoci;
        numHapsInOverlap++;
    }
    window -> startLocus = record -> overlapStartLoci[window -> windowNumOnChromosome % ((record -> WINDOW_SIZE - record -> STEP_SIZE) / record -> STEP_SIZE + 1)];
    window -> endLocus = record -> encoder -> endLocus;
    if (isSameChromosome) {
        record -> numLociInOverlap = window -> numLoci - (record-> HAP_SIZE * record -> STEP_SIZE);
    } else {
        record -> curChrom -> l = 0;
        kputs(ks_str(record -> parser -> nextChromosome), record -> curChrom);
        record -> numLociInOverlap = 0;
        record -> curWinOnChromNum = 0;
    }

    return window;

}

void* sliding_window_multi_threaded(void* arg) {

    WindowRecord* record = (WindowRecord*) arg;

    Window* window;

    double** leftHaps = create_matrix(record -> WINDOW_SIZE, record -> encoder -> numSamples);
    double** rightHaps = create_matrix(record -> WINDOW_SIZE, record -> encoder -> numSamples);

    while (true) {

        pthread_mutex_lock(&genomeLock);
        if (record -> parser -> isEOF) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }

        record -> leftHaps = leftHaps;
        record -> rightHaps = rightHaps;

        window = get_next_window(record);
        
        pthread_mutex_unlock(&genomeLock);

        // Process window.

        pthread_mutex_lock(&listLock);
        *kl_pushp(WindowPtr, record -> windowList) = window;
        pthread_mutex_unlock(&listLock);

    }
    
    destroy_matrix(leftHaps, record -> WINDOW_SIZE);
    destroy_matrix(rightHaps, record -> WINDOW_SIZE);
    
    return NULL;
}

void sliding_window_single_thread(WindowRecord* record) {

    Window* window;
    double** swap;

    while (true) {

        if (record -> parser -> isEOF) {
            break;
        }

        swap = record -> windowCalculations;
        record -> windowCalculations = record -> overlapCalculations;
        record -> overlapCalculations = swap;

        window = get_next_window(record);

        // Process window.

        *kl_pushp(WindowPtr, record -> windowList) = window;

    }

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

    WindowRecord* record = (WindowRecord*) malloc(sizeof(WindowRecord));
    record -> parser = parser;
    record -> encoder = encoder;
    record -> windowList = kl_init(WindowPtr);

    record -> HAP_SIZE = HAP_SIZE;
    record -> STEP_SIZE = STEP_SIZE;
    record -> WINDOW_SIZE = WINDOW_SIZE;
    record -> NUM_THREADS = NUM_THREADS;

    record -> numLociInOverlap = 0;
    record -> overlapStartLoci = (int*) calloc((WINDOW_SIZE - STEP_SIZE) / STEP_SIZE + 1, sizeof(int));
    record -> globalCalculations = create_matrix(encoder -> numSamples, encoder -> numSamples);
    record -> numLociInGlobal = 0;

    record -> curWinNum = 0;
    record -> curWinOnChromNum = 0;
    record -> curChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(ks_str(parser -> nextChromosome), record -> curChrom);

    Window** windows;

    if (NUM_THREADS == 1) {
        record -> overlapLeftHaps = NULL;
        record -> overlapRightHaps = NULL;
        record -> windowCalculations = create_matrix(encoder -> numSamples, encoder -> numSamples);
        record -> overlapCalculations = create_matrix(encoder -> numSamples, encoder -> numSamples);

        sliding_window_single_thread(record);

        destroy_matrix(record -> windowCalculations, encoder -> numSamples);
        destroy_matrix(record -> overlapCalculations, encoder -> numSamples);
    } else {
        record -> windowCalculations = NULL;
        record -> overlapCalculations = NULL;
        record -> leftHaps = NULL;
        record -> rightHaps = NULL;
        record -> overlapLeftHaps = create_matrix(WINDOW_SIZE, encoder -> numSamples);
        record -> overlapRightHaps = create_matrix(WINDOW_SIZE, encoder -> numSamples);

        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, sliding_window_multi_threaded, (void*) record);
    
        sliding_window_multi_threaded((void*) record);

        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);

        free(threads);
        destroy_matrix(record -> overlapLeftHaps, WINDOW_SIZE);
        destroy_matrix(record -> overlapRightHaps, WINDOW_SIZE);
    }

    *numWindows = record -> curWinNum;
    windows = (Window**) malloc(*numWindows * sizeof(Window*));
    int i = 0;
    for (kliter_t(WindowPtr)* it = kl_begin(record -> windowList); it != kl_end(record -> windowList); it = kl_next(it)) {
        windows[i] = kl_val(it);
        i++;
    }
    if (NUM_THREADS > 1)
        sort_windows(windows, *numWindows);

    kl_destroy(WindowPtr, record -> windowList);
    free(record -> overlapStartLoci);
    destroy_matrix(record -> globalCalculations, encoder -> numSamples);
    free(record -> curChrom -> s); free(record -> curChrom);
    free(record);

    return windows;

}