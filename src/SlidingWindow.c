
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

#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

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

    double** winCalcs;
    double** overlapCalcs;
    double** globalCalcs;

    int numLociInOverlap;
    int numLociInGlobal;
    int* overlapStartLoci;

    int curWinNum;
    int curWinOnChromNum;
    kstring_t* curChrom;
} WindowRecord;

Window* get_next_window(WindowRecord* record) {

    Window* window = init_window();
    window -> winNum = ++record -> curWinNum;
    window -> winNumOnChrom = ++record -> curWinOnChromNum;
    kputs(ks_str(record -> curChrom), window -> chromosome);
    window -> numLoci = record -> numLociInOverlap;
    
    int numHapsInOverlap = record -> numLociInOverlap / (record -> HAP_SIZE);
    int numStartsInWin = (record -> WINDOW_SIZE - record -> STEP_SIZE) / record -> STEP_SIZE + 1;
    bool isSameChromosome = true;

    double* swap = NULL;

    /*
    if (record -> NUM_THREADS > 1 && (record -> STEP_SIZE != record -> WINDOW_SIZE) && window -> winNumOnChrom > 1) {
        int overlapSize = record -> WINDOW_SIZE - record -> STEP_SIZE;
        for (int i = 0; i < record -> encoder -> numSamples; i++) {
            for (int j = 0; j < overlapSize; j++) {
                record -> leftHaps[i][j] = record -> overlapLeftHaps[i][j];
                record -> rightHaps[i][j] = record -> overlapRightHaps[i][j];
            }
        }
        for (int i = 0; i < overlapSize; i++) {
            SWAP(record -> overlapLeftHaps[i], record -> overlapLeftHaps[((overlapSize % record -> STEP_SIZE) + i) % overlapSize], swap);
            SWAP(record -> overlapRightHaps[i], record -> overlapRightHaps[((overlapSize % record -> STEP_SIZE) + i) % overlapSize], swap);
        }
    }
    */

    while(numHapsInOverlap < record -> WINDOW_SIZE && isSameChromosome) {
        isSameChromosome = get_next_haplotype(record -> parser, record -> encoder, true, record -> HAP_SIZE);

        if (record -> NUM_THREADS == 1) {
            // TODO.
        } else {
            SWAP(record -> leftHaps[numHapsInOverlap], record -> encoder -> leftHaplotype, swap);
            SWAP(record -> rightHaps[numHapsInOverlap], record -> encoder -> rightHaplotype, swap);
            if (numHapsInOverlap >= record -> STEP_SIZE) {
                for (int i = 0; i < record -> encoder -> numSamples; i++) {
                    record -> overlapLeftHaps[numHapsInOverlap - record -> STEP_SIZE][i] = record -> leftHaps[numHapsInOverlap][i];
                    record -> overlapRightHaps[numHapsInOverlap - record -> STEP_SIZE][i] = record -> rightHaps[numHapsInOverlap][i];
                }
            }
            printf("\n");
            printf("Window %d:\n", record -> winNum);
            for (int i = 0; i < record -> encoder -> numSamples; i++) {
                for (int j = 0; j < record -> WINDOW_SIZE; j++) {
                    printf("%5f/%5f\t", record -> leftHaps[i][j], record -> rightHaps[i][j]);
                }
                printf("\n");
            }
            printf("\nOverlap:\n");
            for (int i = 0; i < record -> encoder -> numSamples; i++) {
                for (int j = 0; j < record -> WINDOW_SIZE - record -> STEP_SIZE; j++) {
                    printf("%5f/%5f\t", record -> overlapLeftHaps[i][j], record -> overlapRightHaps[i][j]);
                }
                printf("\n");
            }
        }

        if (numHapsInOverlap % record -> STEP_SIZE == 0)
            record -> overlapStartLoci[(window -> winNumOnChrom + numHapsInOverlap) % numStartsInWin] = record -> encoder -> startLocus;
        window -> numLoci += record -> encoder -> numLoci;
        numHapsInOverlap++;
    }
    window -> startLocus = record -> overlapStartLoci[(window -> winNumOnChrom + numHapsInOverlap) % numStartsInWin];
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

        SWAP(record -> winCalcs, record -> overlapCalcs, swap);

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
        while (j >= 0 && windows[j] -> winNum > temp -> winNum) {
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
    record -> globalCalcs = create_matrix(encoder -> numSamples, encoder -> numSamples);
    record -> numLociInGlobal = 0;

    record -> curWinNum = 0;
    record -> curWinOnChromNum = 0;
    record -> curChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(ks_str(parser -> nextChromosome), record -> curChrom);

    Window** windows;

    if (NUM_THREADS == 1) {
        record -> overlapLeftHaps = NULL;
        record -> overlapRightHaps = NULL;
        record -> winCalcs = create_matrix(encoder -> numSamples, encoder -> numSamples);
        record -> overlapCalcs = create_matrix(encoder -> numSamples, encoder -> numSamples);

        sliding_window_single_thread(record);

        destroy_matrix(record -> winCalcs, encoder -> numSamples);
        destroy_matrix(record -> overlapCalcs, encoder -> numSamples);
    } else {
        record -> winCalcs = NULL;
        record -> overlapCalcs = NULL;
        record -> leftHaps = NULL;
        record -> rightHaps = NULL;
        record -> overlapLeftHaps = create_matrix(WINDOW_SIZE - STEP_SIZE, encoder -> numSamples);
        record -> overlapRightHaps = create_matrix(WINDOW_SIZE - STEP_SIZE, encoder -> numSamples);

        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, sliding_window_multi_threaded, (void*) record);
    
        sliding_window_multi_threaded((void*) record);

        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);

        free(threads);
        destroy_matrix(record -> overlapLeftHaps, WINDOW_SIZE - STEP_SIZE);
        destroy_matrix(record -> overlapRightHaps, WINDOW_SIZE - STEP_SIZE);
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
    destroy_matrix(record -> globalCalcs, encoder -> numSamples);
    free(record -> curChrom -> s); free(record -> curChrom);
    free(record);

    return windows;

}