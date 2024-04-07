
// File: SlidingWindow.c
// Date: 
// Author: TQ Smith
// Purpose: 

#include "SlidingWindow.h"

#include "../lib/klist.h"
#define do_not_free_window_ptr(w)
KLIST_INIT(WindowPtr, Window*, do_not_free_window_ptr)

#include <pthread.h>

#include "AlleleSharingDistance.h"

#include "Matrix.h"
MATRIX_INIT(double, double)
MATRIX_INIT(ibs, IBS)
MATRIX_INIT(geno, Genotype)

#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t windowListLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    VCFLocusParser* parser;
    HaplotypeEncoder* encoder;
    klist_t(WindowPtr)* winList;

    int HAP_SIZE;
    int STEP_SIZE;
    int WINDOW_SIZE;
    int numSamples;

    Genotype** winGeno;
    int numHapsInOverlap;
    int winStartIndex;
    int winEndIndex;
    int* winStartLoci;

    int globalNumHaps;

    int curWinNum;
    int curWinOnChromNum;
    kstring_t* curChrom;
} WindowRecord;

Window* get_next_window(WindowRecord* record, Genotype** threadGeno) {

    Window* window = init_window();
    window -> winNum = record -> curWinNum++;
    window -> winNumOnChrom = record -> curWinOnChromNum++;
    kputs(ks_str(record -> curChrom), window -> chromosome);
    window -> numLoci = record -> numHapsInOverlap * record -> HAP_SIZE;

    bool isSameChrom = true;
    Genotype* temp;

    if (window -> winNumOnChrom != 1)
        record -> winStartIndex = (record -> winStartIndex + record -> STEP_SIZE) % record -> WINDOW_SIZE;

    while(record -> numHapsInOverlap < record -> WINDOW_SIZE && isSameChrom) {
        isSameChrom = get_next_haplotype(record -> parser, record -> encoder, true, record -> HAP_SIZE);
        if (threadGeno == NULL) {
            SWAP(record -> winGeno[record -> winEndIndex], record -> encoder -> genotypes, temp);
        } else {
            SWAP(record -> winGeno[record -> numHapsInOverlap], record -> encoder -> genotypes, temp);
            if (isSameChrom && record -> numHapsInOverlap >= record -> STEP_SIZE)
                for (int i = 0; i < record -> numSamples; i++)
                    threadGeno[record -> numHapsInOverlap][i] = record -> winGeno[record -> winEndIndex][i];
        }
        record -> winEndIndex = (record -> winEndIndex + 1) % record -> WINDOW_SIZE;
        if (record -> numHapsInOverlap % record -> STEP_SIZE == 0)
            record -> winStartLoci[(window -> winNumOnChrom + record -> numHapsInOverlap / record -> STEP_SIZE) % (record -> WINDOW_SIZE / record -> STEP_SIZE + 1)] = record -> encoder -> startLocus;
        window -> numLoci += record -> encoder -> numLoci;
        record -> numHapsInOverlap++;
    }

    window -> startLocus = record -> winStartLoci[window -> winNumOnChrom % (record -> WINDOW_SIZE / record -> STEP_SIZE + 1)];
    window -> endLocus = record -> encoder -> endLocus;
    if (isSameChrom) {
        record -> numHapsInOverlap = record -> WINDOW_SIZE - record -> STEP_SIZE;
    } else {
        record -> curChrom -> l = 0;
        kputs(ks_str(record -> parser -> nextChrom), record -> curChrom);
        record -> numHapsInOverlap = 0;
        record -> curWinOnChromNum = 1;
        record -> winStartIndex = 0;
        record -> winEndIndex = 0;
    }

    return window;

}

void* sliding_window_multi_thread(void* arg) {
    return NULL;
}

void sliding_window_single_thread(WindowRecord* record) {

    Window* window;

    while (!(record -> parser -> isEOF)) {

        window = get_next_window(record, NULL);

        if (record -> curWinOnChromNum == 1)
            record -> globalNumHaps += window -> numLoci / record -> HAP_SIZE;
        else
            record -> globalNumHaps += record -> STEP_SIZE;

        printf("Window %d Contents:\n", window -> winNum);
        int start = record -> winStartIndex;
        for (int i = 0; i < record -> WINDOW_SIZE; i++) {
            for (int j = 0; j < record -> numSamples; j++) {
                printf("%3ld/%-3ld\t", record -> winGeno[start][j].left, record -> winGeno[start][j].right);
            }
            printf("\n");
            start = (start + 1) % record -> WINDOW_SIZE;
        }
        printf("\n");

        *kl_pushp(WindowPtr, record -> winList) = window;

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

Window** sliding_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows) {
    
    WindowRecord* record = (WindowRecord*) calloc(1, sizeof(WindowRecord));
    record -> parser = parser;
    record -> encoder = encoder;
    record -> winList = kl_init(WindowPtr);

    record -> HAP_SIZE = HAP_SIZE;
    record -> STEP_SIZE = STEP_SIZE;
    record -> WINDOW_SIZE = WINDOW_SIZE;
    record -> numSamples = encoder -> numSamples;

    record -> winGeno = create_matrix(geno, WINDOW_SIZE, encoder -> numSamples);
    record -> numHapsInOverlap = 0;
    record -> winStartIndex = 0;
    record -> winEndIndex = 0;
    record -> winStartLoci = (int*) calloc(WINDOW_SIZE / STEP_SIZE + 1, sizeof(int));

    record -> globalNumHaps = 0;

    record -> curWinNum = 1;
    record -> curWinOnChromNum = 1;
    record -> curChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(ks_str(parser -> nextChrom), record -> curChrom);

    Window** windows;

    if (NUM_THREADS == 1) {
        sliding_window_single_thread(record);
    } else {
        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, sliding_window_multi_thread, (void*) record);
    
        sliding_window_multi_thread((void*) record);

        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);
    }
    
    Window* global = init_window();
    kputs("Global", global -> chromosome);
    global -> winNum = 0;
    global -> winNumOnChrom = 0;
    global -> startLocus = 0;
    global -> endLocus = 0;
    global -> numLoci = record -> globalNumHaps;
    
    *numWindows = record -> curWinNum;
    windows = (Window**) malloc(*numWindows * sizeof(Window*));
    int i = 1;
    for (kliter_t(WindowPtr)* it = kl_begin(record -> winList); it != kl_end(record -> winList); it = kl_next(it)) {
        windows[i] = kl_val(it);
        i++;
    }
    windows[0] = global;
    if (NUM_THREADS > 1)
        sort_windows(windows, *numWindows);

    kl_destroy(WindowPtr, record -> winList);
    destroy_matrix(geno, record -> winGeno, WINDOW_SIZE);
    free(record -> winStartLoci);
    free(record -> curChrom -> s); free(record -> curChrom);
    free(record);

    return windows;

}

typedef struct {
    VCFLocusParser* parser;
    HaplotypeEncoder* encoder;
    IBS** alleleCounts;
    int numLoci;
    int HAP_SIZE;
} GlobalWindowRecord;

void* global_window_multi_thread(void* arg) {
    GlobalWindowRecord* record = (GlobalWindowRecord*) arg;
    int numSamples = record -> encoder -> numSamples;

    IBS** alleleCounts = create_matrix(ibs, numSamples, numSamples);
    Genotype* genotypes = (Genotype*) calloc(numSamples, sizeof(Genotype));
    Genotype* temp;

    pthread_mutex_lock(&genomeLock);
    while (!record -> parser -> isEOF) {
        get_next_haplotype(record -> parser, record -> encoder, true, record -> HAP_SIZE);
        record -> numLoci += record -> encoder -> numLoci;
        SWAP(genotypes, record -> encoder -> genotypes, temp);
        pthread_mutex_unlock(&genomeLock);

        pairwise_ibs(alleleCounts, genotypes, numSamples);
    }
    pthread_mutex_unlock(&genomeLock);

    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < numSamples; i++) {
        for (int j = i + 1; j < numSamples; j++) {
            record -> alleleCounts[i][j].ibs0 += alleleCounts[i][j].ibs0;
            record -> alleleCounts[i][j].ibs1 += alleleCounts[i][j].ibs1;
            record -> alleleCounts[i][j].ibs2 += alleleCounts[i][j].ibs2;
        }
    }
    pthread_mutex_unlock(&globalLock);

    destroy_matrix(ibs, alleleCounts, numSamples);
    free(genotypes);
    return NULL;
}

Window* global_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int NUM_THREADS) {

    IBS** alleleCounts = create_matrix(ibs, encoder -> numSamples, encoder -> numSamples);
    double** asd = create_matrix(double, encoder -> numSamples, encoder -> numSamples);

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
            pairwise_ibs(alleleCounts, encoder -> genotypes, encoder -> numSamples);
        }
    } else {
        GlobalWindowRecord* record = (GlobalWindowRecord*) calloc(1, sizeof(GlobalWindowRecord));
        record -> parser = parser;
        record -> encoder = encoder;
        record -> alleleCounts = alleleCounts;
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

    for (int i = 0; i < encoder -> numSamples; i++) {
        for (int j = i + 1; j < encoder -> numSamples; j++) {
            asd[i][j] = asd[j][i] = ibs_to_asd(alleleCounts[i][j]);
        }
    }

    destroy_matrix(ibs, alleleCounts, encoder -> numSamples);
    destroy_matrix(double, asd, encoder -> numSamples);

    return window;

}