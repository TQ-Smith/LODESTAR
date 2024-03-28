
// File: SlidingWindow.c
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Slides a window along the contents of a VCF file.

#include "SlidingWindow.h"

#include "Matrix.h"

#include "AlleleSharingDistance.h"

#include "../klib/klist.h"
#define do_not_free_window_ptr(w)
KLIST_INIT(WindowPtr, Window*, do_not_free_window_ptr)

#include <pthread.h>

#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t listLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    VCFGenotypeParser* parser;
    HaplotypeEncoder* encoder;
    klist_t(WindowPtr)* windowList;

    int HAP_SIZE;
    int STEP_SIZE;
    int WINDOW_SIZE;
    int NUM_THREADS;
    int numSamples;

    double** leftHaps;
    double** rightHaps;
    double** overlapLeftHaps;
    double** overlapRightHaps;
    double** globalIBS;
    int globalNumLoci;
    int startIndex;
    int endIndex;
    int numLociInOverlap;
    int* overlapStartLoci;
    bool isSameChrom;

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
    
    int numHapsInWin = record -> numLociInOverlap / (record -> HAP_SIZE);
    int numStartsInWin = (record -> WINDOW_SIZE - record -> STEP_SIZE) / record -> STEP_SIZE + 2;
    int overlapSize = record -> WINDOW_SIZE - record -> STEP_SIZE;
    bool isSameChrom = true;

    double* swap = NULL;

    if (window -> winNumOnChrom == 1) {
        record -> startIndex = 0;
        record -> endIndex = 0;
    } else if (record -> NUM_THREADS == 1) {
        record -> startIndex = (record -> startIndex + record -> STEP_SIZE) % record -> WINDOW_SIZE;
    } else {
        int start = record -> startIndex;
        for (int i = 0; i < overlapSize; i++) {
            for (int j = 0; j < record -> numSamples; j++) {
                record -> leftHaps[i][j] = record -> overlapLeftHaps[start][j];
                record -> rightHaps[i][j] = record -> overlapRightHaps[start][j];
            }
            start = (start + 1) % overlapSize;
        }
        record -> startIndex = (record -> startIndex + record -> STEP_SIZE) % overlapSize;
    }

    while(numHapsInWin < record -> WINDOW_SIZE && isSameChrom) {
        isSameChrom = get_next_haplotype(record -> parser, record -> encoder, true, record -> HAP_SIZE);

        if (record -> NUM_THREADS == 1) {
            SWAP(record -> leftHaps[record -> endIndex], record -> encoder -> leftHaps, swap);
            SWAP(record -> rightHaps[record -> endIndex], record -> encoder -> rightHaps, swap);
            record -> endIndex = (record -> endIndex + 1) % record -> WINDOW_SIZE;
        } else {
            SWAP(record -> leftHaps[numHapsInWin], record -> encoder -> leftHaps, swap);
            SWAP(record -> rightHaps[numHapsInWin], record -> encoder -> rightHaps, swap);
            if (isSameChrom && numHapsInWin >= record -> STEP_SIZE) {
                for (int i = 0; i < record -> numSamples; i++) {
                    record -> overlapLeftHaps[record -> endIndex][i] = record -> leftHaps[numHapsInWin][i];
                    record -> overlapRightHaps[record -> endIndex][i] = record -> rightHaps[numHapsInWin][i];
                } 
                record -> endIndex = (record -> endIndex + 1) % overlapSize;
            }
        }

        if (numHapsInWin % record -> STEP_SIZE == 0)
            record -> overlapStartLoci[(window -> winNumOnChrom + numHapsInWin / record -> STEP_SIZE) % numStartsInWin] = record -> encoder -> startLocus;
        window -> numLoci += record -> encoder -> numLoci;
        numHapsInWin++;
    }
    window -> startLocus = record -> overlapStartLoci[window -> winNumOnChrom % numStartsInWin];
    window -> endLocus = record -> encoder -> endLocus;
    record -> isSameChrom = isSameChrom;
    if (isSameChrom) {
        record -> numLociInOverlap = window -> numLoci - (record-> HAP_SIZE * record -> STEP_SIZE);
    } else {
        record -> curChrom -> l = 0;
        kputs(ks_str(record -> parser -> nextChrom), record -> curChrom);
        record -> numLociInOverlap = 0;
        record -> curWinOnChromNum = 0;
    }

    return window;

}

void* sliding_window_multi_thread(void* arg) {

    WindowRecord* record = (WindowRecord*) arg;

    Window* window;

    double** leftHaps = create_matrix(record -> WINDOW_SIZE, record -> numSamples);
    double** rightHaps = create_matrix(record -> WINDOW_SIZE, record -> numSamples);
    double** winIBS = create_matrix(record -> numSamples, record -> numSamples);
    double** globalIBS = create_matrix(record -> numSamples, record -> numSamples);

    int globalNumLoci = 0;

    bool isSameChrom;

    while (true) {

        pthread_mutex_lock(&genomeLock);
        if (record -> parser -> isEOF) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        record -> leftHaps = leftHaps;
        record -> rightHaps = rightHaps;

        window = get_next_window(record);
        isSameChrom = record -> isSameChrom;
        pthread_mutex_unlock(&genomeLock);

        int numHapsInWin = (int) ceil((double) window -> numLoci / record -> HAP_SIZE);
        
        for (int i = 0; i < numHapsInWin; i++) {
            process_haplotype_multi_thread(leftHaps[i], rightHaps[i], winIBS, globalIBS, record -> numSamples, record -> STEP_SIZE, numHapsInWin, isSameChrom, i);
        }

        if (isSameChrom) {
            globalNumLoci += record -> STEP_SIZE;
        } else {
            globalNumLoci += numHapsInWin;
        }

        printf("Window %d Contents:\n", window -> winNum);
        for (int i = 0; i < numHapsInWin; i++) {
            for (int j = 0; j < record -> numSamples; j++) {
                printf("%5f/%5f\t", leftHaps[i][j], rightHaps[i][j]);
            }
            printf("\n");
        }
        printf("Window %d ASD:\n", window -> winNum);
        for (int i = 0; i < record -> numSamples; i++) {
            for (int j = 0; j < record -> numSamples; j++) {
                printf("%5f\t", winIBS[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        pthread_mutex_lock(&listLock);
        *kl_pushp(WindowPtr, record -> windowList) = window;
        pthread_mutex_unlock(&listLock);

    }

    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < record -> numSamples; i++) {
        for (int j = i + 1; j < record -> numSamples; j++) {
            record -> globalIBS[i][j] += globalIBS[i][j];
            record -> globalIBS[j][i] += globalIBS[j][i];
        }
    }
    record -> globalNumLoci += globalNumLoci;
    pthread_mutex_unlock(&globalLock);
    
    destroy_matrix(leftHaps, record -> WINDOW_SIZE);
    destroy_matrix(rightHaps, record -> WINDOW_SIZE);
    destroy_matrix(winIBS, record -> numSamples);
    destroy_matrix(globalIBS, record -> numSamples);
    
    return NULL;
}

void sliding_window_single_thread(WindowRecord* record) {

    Window* window;

    double** leftHaps = create_matrix(record -> WINDOW_SIZE, record -> numSamples);
    double** rightHaps = create_matrix(record -> WINDOW_SIZE, record -> numSamples);
    double** winIBS = create_matrix(record -> numSamples, record -> numSamples);
    double** offsetIBS = create_matrix(record -> numSamples, record -> numSamples);
    double** winASD = create_matrix(record -> numSamples, record -> numSamples);
    record -> leftHaps = leftHaps;
    record -> rightHaps = rightHaps;

    int start;

    while (true) {

        if (record -> parser -> isEOF) {
            break;
        }

        window = get_next_window(record);

        int numHapsInWin = (int) ceil((double) window -> numLoci / record -> HAP_SIZE);

        if (window -> winNumOnChrom == 1) {
            start = record -> startIndex;
            for (int i = 0; i < numHapsInWin; i++) {
                process_haplotype_single_thread(leftHaps[start], rightHaps[start], winIBS, offsetIBS, winASD, record -> globalIBS, record -> numSamples, record -> STEP_SIZE, numHapsInWin, record -> isSameChrom, i);
                start = (start + 1) % record -> WINDOW_SIZE;
            }
        } else {

        }

        if (record -> isSameChrom) {
            record -> globalNumLoci += record -> STEP_SIZE;
        } else {
            record -> globalNumLoci += numHapsInWin;
        }

        printf("Window %d Contents:\n", window -> winNum);
        start = record -> startIndex;
        for (int i = 0; i < numHapsInWin; i++) {
            for (int j = 0; j < record -> numSamples; j++) {
                printf("%5f/%5f\t", leftHaps[start][j], rightHaps[start][j]);
            }
            printf("\n");
            start = (start + 1) % record -> WINDOW_SIZE;
        }
        printf("Window %d ASD:\n", window -> winNum);
        for (int i = 0; i < record -> numSamples; i++) {
            for (int j = 0; j < record -> numSamples; j++) {
                printf("%5f\t", winIBS[i][j]);
            }
            printf("\n");
        }
        printf("\n");

        *kl_pushp(WindowPtr, record -> windowList) = window;

    }

    destroy_matrix(leftHaps, record -> WINDOW_SIZE);
    destroy_matrix(rightHaps, record -> WINDOW_SIZE);
    destroy_matrix(winIBS, record -> numSamples);
    destroy_matrix(offsetIBS, record -> numSamples);
    destroy_matrix(winASD, record -> numSamples);

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
    record -> numSamples = encoder -> numSamples;

    record -> overlapStartLoci = (int*) calloc((WINDOW_SIZE - STEP_SIZE) / STEP_SIZE + 2, sizeof(int));
    record -> globalIBS = create_matrix(record -> numSamples, record -> numSamples);
    record -> globalNumLoci = 0;

    record -> numLociInOverlap = 0;
    record -> startIndex = 0;
    record -> endIndex = 0;

    record -> curWinNum = 0;
    record -> curWinOnChromNum = 0;
    record -> curChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(ks_str(parser -> nextChrom), record -> curChrom);

    Window** windows;

    if (NUM_THREADS == 1) {
        sliding_window_single_thread(record);
    } else {
        record -> overlapLeftHaps = create_matrix(WINDOW_SIZE - STEP_SIZE, encoder -> numSamples);
        record -> overlapRightHaps = create_matrix(WINDOW_SIZE - STEP_SIZE, encoder -> numSamples);
        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, sliding_window_multi_thread, (void*) record);
    
        sliding_window_multi_thread((void*) record);

        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);
        destroy_matrix(record -> overlapLeftHaps, WINDOW_SIZE - STEP_SIZE);
        destroy_matrix(record -> overlapRightHaps, WINDOW_SIZE - STEP_SIZE);
        free(threads);
    }

    Window* global = init_window();
    global -> winNum = 0;
    global -> winNumOnChrom = 0;
    global -> startLocus = 0;
    global -> endLocus = 0;
    kputs("Global", global -> chromosome);
    global -> numLoci = record -> globalNumLoci;
    for (int i = 0; i < record -> numSamples; i++) {
        for (int j = i + 1; j < record -> numSamples; j++) {
            record -> globalIBS[i][j] = record -> globalIBS[j][i] = 1.0 - (record -> globalIBS[i][j] / (2 * record -> globalIBS[j][i]));
        }
    }

    printf("Global:\n");
    for (int i = 0; i < record -> numSamples; i++) {
        for (int j = 0; j < record -> numSamples; j++) {
            printf("%5f\t", record -> globalIBS[i][j]);
        }
        printf("\n");
    }
    
    *numWindows = record -> curWinNum;
    windows = (Window**) malloc((*numWindows + 1) * sizeof(Window*));
    int i = 0;
    for (kliter_t(WindowPtr)* it = kl_begin(record -> windowList); it != kl_end(record -> windowList); it = kl_next(it)) {
        windows[i] = kl_val(it);
        i++;
    }
    windows[*numWindows] = global;
    *numWindows = *numWindows + 1;
    if (NUM_THREADS > 1)
        sort_windows(windows, *numWindows);
    
    kl_destroy(WindowPtr, record -> windowList);
    free(record -> overlapStartLoci);
    destroy_matrix(record -> globalIBS, record -> numSamples);
    free(record -> curChrom -> s); free(record -> curChrom);
    free(record);

    return windows;

}