
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
MATRIX_INIT(geno, Genotype)
MATRIX_INIT(double, double)

#include <math.h>

#include "MultidimensionalScaling.h"

#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t windowListLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

typedef struct {
    VCFLocusParser* parser;
    HaplotypeEncoder* encoder;
    RealSymEigen* eigen;
    klist_t(WindowPtr)* winList;
    int k;

    int HAP_SIZE;
    int STEP_SIZE;
    int WINDOW_SIZE;
    int numSamples;

    Genotype** winGeno;
    int numHapsInOverlap;
    int winStartIndex;
    int winEndIndex;
    int* winStartLoci;

    IBS* globalAlleleCounts;
    int globalNumHaps;

    int curWinNum;
    int curWinNumOnChrom;
    kstring_t* curChrom;
} WindowRecord;

void perform_mds_on_window(Window* window, RealSymEigen* eigen, double* asd, int k) {
    double** X = create_matrix(double, eigen -> N, k);
    int INFO = compute_classical_mds(eigen, asd, k, X);
    if (INFO != 0) {
        destroy_matrix(double, X, eigen -> N);
        return;
    }
    window -> X = X;
}

Window* get_next_window(WindowRecord* record, Genotype** threadGeno) {

    Window* window = init_window();
    window -> winNum = record -> curWinNum++;
    window -> winNumOnChrom = record -> curWinNumOnChrom++;
    ks_overwrite(ks_str(record -> curChrom), window -> chromosome);
    window -> numLoci = record -> numHapsInOverlap * record -> HAP_SIZE;

    bool isSameChrom = true;
    Genotype* temp;

    if (window -> winNumOnChrom == 1) {
        record -> winStartIndex = 0;
        record -> winEndIndex = 0;
    } else {
        record -> winStartIndex = (record -> winStartIndex + record -> STEP_SIZE) % record -> WINDOW_SIZE;
    }

    while(record -> numHapsInOverlap < record -> WINDOW_SIZE && isSameChrom) {
        isSameChrom = get_next_haplotype(record -> parser, record -> encoder, record -> HAP_SIZE);
        if (threadGeno == NULL) {
            SWAP(record -> winGeno[record -> winEndIndex], record -> encoder -> genotypes, temp);
        } else {
            SWAP(record -> winGeno[record -> winEndIndex], record -> encoder -> genotypes, temp);
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
        record -> globalNumHaps += record -> STEP_SIZE;
        record -> numHapsInOverlap = record -> WINDOW_SIZE - record -> STEP_SIZE;
    } else {
        record -> globalNumHaps += window -> numLoci / record -> HAP_SIZE;
        record -> curChrom -> l = 0;
        ks_overwrite(ks_str(record -> parser -> nextChrom), record -> curChrom);
        record -> numHapsInOverlap = 0;
        record -> curWinNumOnChrom = 1;
    }

    return window;

}

void* sliding_window_multi_thread(void* arg) {

    WindowRecord* record = (WindowRecord*) arg;
    int numSamples = record -> numSamples;
    int HAP_SIZE = record -> HAP_SIZE;
    int STEP_SIZE = record -> STEP_SIZE;
    int k = record -> k;

    Window* window;

    Genotype** threadGeno = create_matrix(geno, record -> WINDOW_SIZE, record -> numSamples);
    IBS* winAlleleCounts = (IBS*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS));
    IBS* globalAlleleCounts = (IBS*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS));
    double* asd = (double*) calloc(PACKED_SIZE(record -> numSamples), sizeof(double));
    RealSymEigen* eigen = init_real_sym_eigen(record -> numSamples);

    int start, numHapsInWin;
    bool isLastWinOnChrom;

    while (true) {

        pthread_mutex_lock(&genomeLock);
        if (record -> parser -> isEOF) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }

        if (record -> curWinNumOnChrom != 1) {
            start = (record -> winStartIndex + record -> STEP_SIZE) % record -> WINDOW_SIZE;
            for (int i = 0; i < record -> WINDOW_SIZE - record -> STEP_SIZE; i++) {
                for (int j = 0; j < record -> numSamples; j++)
                    threadGeno[i][j] = record -> winGeno[start][j];
                start = (start + 1) % record -> WINDOW_SIZE;
            }
        }

        window = get_next_window(record, threadGeno);

        isLastWinOnChrom = record -> curWinNumOnChrom == 1;

        pthread_mutex_unlock(&genomeLock);

        numHapsInWin = (int) ceil((double) window -> numLoci / HAP_SIZE);

        process_window_multi_thread(threadGeno, winAlleleCounts, globalAlleleCounts, asd, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE);

        printf("Window %d ASD:\n", window -> winNum);
        for (int i = 0; i < record -> numSamples; i++) {
            for (int j = i + 1; j < record -> numSamples; j++) {
                printf("%5f\t", asd[PACKED_INDEX(i, j)]);
            }
            printf("\n");
        }
        printf("\n");

        perform_mds_on_window(window, eigen, asd, k);

        pthread_mutex_lock(&windowListLock);
        *kl_pushp(WindowPtr, record -> winList) = window;
        pthread_mutex_unlock(&windowListLock);

    }

    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < record -> numSamples; i++)
        for (int j = i + 1; j < record -> numSamples; j++) 
            add_ibs(&(record -> globalAlleleCounts[PACKED_INDEX(i, j)]), &globalAlleleCounts[PACKED_INDEX(i, j)]);
    pthread_mutex_unlock(&globalLock);

    destroy_matrix(geno, threadGeno, record -> WINDOW_SIZE);
    free(winAlleleCounts);
    free(globalAlleleCounts);
    free(asd);
    destroy_real_sym_eigen(eigen);

    return NULL;

}

void sliding_window_single_thread(WindowRecord* record) {

    Window* window;

    IBS* winAlleleCounts = (IBS*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS));
    IBS* stepAlleleCounts = (IBS*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS));
    double* asd = (double*) calloc(PACKED_SIZE(record -> numSamples), sizeof(double));

    int numHapsInWin;

    while (!(record -> parser -> isEOF)) {

        window = get_next_window(record, NULL);
        
        numHapsInWin = (int) ceil((double) window -> numLoci / record -> HAP_SIZE);

        /*
        printf("Window %d:\n", window -> winNum);
        int start = record -> winStartIndex;
        for(int i = 0; i < record -> WINDOW_SIZE; i++) {
            for (int j = 0; j < record -> numSamples; j++)
                printf("%lu/%lu\t", record -> winGeno[start][j].left, record -> winGeno[start][j].right);
            start = (start + 1) % record -> WINDOW_SIZE;
            printf("\n");
        }
        */
        
        process_window_single_thread(record -> winGeno, record -> winStartIndex, winAlleleCounts, stepAlleleCounts, record -> globalAlleleCounts, asd, numHapsInWin, window -> winNumOnChrom == 1, record -> curWinNumOnChrom == 1, record -> numSamples, record -> STEP_SIZE, record -> WINDOW_SIZE);
        
        printf("Window %d ASD:\n", window -> winNum);
        for (int i = 0; i < record -> numSamples; i++) {
            for (int j = i + 1; j < record -> numSamples; j++) {
                printf("%5f\t", asd[PACKED_INDEX(i, j)]);
            }
            printf("\n");
        }
        printf("\n");

        perform_mds_on_window(window, record -> eigen, asd, record -> k);

        *kl_pushp(WindowPtr, record -> winList) = window;

    }

    free(winAlleleCounts);
    free(stepAlleleCounts);
    free(asd);

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

Window** sliding_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int k, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows) {
    
    WindowRecord* record = (WindowRecord*) calloc(1, sizeof(WindowRecord));
    record -> parser = parser;
    record -> encoder = encoder;
    record -> eigen = init_real_sym_eigen(encoder -> numSamples);
    record -> winList = kl_init(WindowPtr);

    record -> k = k;

    record -> HAP_SIZE = HAP_SIZE;
    record -> STEP_SIZE = STEP_SIZE;
    record -> WINDOW_SIZE = WINDOW_SIZE;
    record -> numSamples = encoder -> numSamples;

    record -> winGeno = create_matrix(geno, WINDOW_SIZE, encoder -> numSamples);
    record -> numHapsInOverlap = 0;
    record -> winStartIndex = 0;
    record -> winEndIndex = 0;
    record -> winStartLoci = (int*) calloc(WINDOW_SIZE / STEP_SIZE + 1, sizeof(int));

    record -> globalAlleleCounts = (IBS*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS));
    record -> globalNumHaps = 0;

    record -> curWinNum = 1;
    record -> curWinNumOnChrom = 1;
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
    ks_overwrite("Global", global -> chromosome);
    global -> winNum = 0;
    global -> winNumOnChrom = 0;
    global -> startLocus = 0;
    global -> endLocus = 0;
    global -> numLoci = record -> globalNumHaps;

    printf("Global ASD:\n");
    for (int i = 0; i < record -> numSamples; i++) {
        for (int j = i + 1; j < record -> numSamples; j++) {
            printf("%5f\t", ibs_to_asd(record -> globalAlleleCounts[PACKED_INDEX(i, j)]));
        }
        printf("\n");
    }

    double* globalASD = (double*) malloc(PACKED_SIZE(record -> numSamples) * sizeof(double));
    for (int i = 0; i < record -> numSamples; i++) 
        for (int j = i + 1; j < record -> numSamples; j++)
            globalASD[PACKED_INDEX(i, j)] = ibs_to_asd(record -> globalAlleleCounts[PACKED_INDEX(i, j)]);
    perform_mds_on_window(global, record -> eigen, globalASD, k);
    free(globalASD);

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
    destroy_real_sym_eigen(record -> eigen);
    free(record -> winStartLoci);
    free(record -> globalAlleleCounts);
    free(record -> curChrom -> s); free(record -> curChrom);
    free(record);

    return windows;

}

typedef struct {
    VCFLocusParser* parser;
    HaplotypeEncoder* encoder;
    IBS* alleleCounts;
    int numLoci;
    int HAP_SIZE;
} GlobalWindowRecord;

void* global_window_multi_thread(void* arg) {
    GlobalWindowRecord* record = (GlobalWindowRecord*) arg;
    int numSamples = record -> encoder -> numSamples;

    IBS* alleleCounts = (IBS*) calloc(PACKED_SIZE(numSamples), sizeof(IBS));
    Genotype* genotypes = (Genotype*) calloc(numSamples, sizeof(Genotype));
    Genotype* temp;

    pthread_mutex_lock(&genomeLock);
    while (!record -> parser -> isEOF) {
        get_next_haplotype(record -> parser, record -> encoder, record -> HAP_SIZE);
        record -> numLoci += record -> encoder -> numLoci;
        SWAP(genotypes, record -> encoder -> genotypes, temp);
        pthread_mutex_unlock(&genomeLock);

        pairwise_ibs(genotypes, alleleCounts, numSamples);
    }
    pthread_mutex_unlock(&genomeLock);

    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < numSamples; i++)
        for (int j = i + 1; j < numSamples; j++)
            add_ibs(&(record -> alleleCounts[PACKED_INDEX(i, j)]), &alleleCounts[PACKED_INDEX(i, j)]);
    pthread_mutex_unlock(&globalLock);

    free(alleleCounts);
    free(genotypes);
    return NULL;
}

Window* global_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int k, int HAP_SIZE, int NUM_THREADS) {

    IBS* alleleCounts = (IBS*) calloc(PACKED_SIZE(encoder -> numSamples), sizeof(IBS));
    double* asd = (double*) calloc(PACKED_SIZE(encoder -> numSamples), sizeof(double));

    Window* window = init_window();
    window -> winNum = 0;
    window -> winNumOnChrom = 0;
    ks_overwrite("Global", window -> chromosome);
    window -> startLocus = 0;
    window -> endLocus = 0;
    window -> asd = (double*) calloc(PACKED_SIZE(encoder -> numSamples), sizeof(double));
    window -> saveIBS = true;
    window -> saveASD = true;

    if (NUM_THREADS == 1) {
        while (!parser -> isEOF) {
            get_next_haplotype(parser, encoder, HAP_SIZE);
            window -> numLoci += encoder -> numLoci;
            pairwise_ibs(encoder -> genotypes, alleleCounts, encoder -> numSamples);
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
            asd[PACKED_INDEX(i, j)] = ibs_to_asd(alleleCounts[PACKED_INDEX(i, j)]);
            window -> asd[PACKED_INDEX(i, j)] = asd[PACKED_INDEX(i, j)];
        }
    }

    RealSymEigen* eigen = init_real_sym_eigen(encoder -> numSamples);
    perform_mds_on_window(window, eigen, asd, k);
    destroy_real_sym_eigen(eigen);

    window -> ibs = alleleCounts;

    free(asd);

    return window;

}