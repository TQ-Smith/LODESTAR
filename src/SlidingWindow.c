
// File: SlidingWindow.c
// Date: 8 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Perform ASD and MDS along windows in the genome or just genome-wide.

#include "SlidingWindow.h"

#include "klist.h"
#define do_not_free_window_ptr(w)
KLIST_INIT(WindowPtr, Window_t*, do_not_free_window_ptr)

#include "Matrix.h"
MATRIX_INIT(geno, Genotype_t)
MATRIX_INIT(double, double)

#include "Logger.h"
#include "AlleleSharingDistance.h"
#include "MultidimensionalScaling.h"
#include <pthread.h>
#include <math.h>

#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

// Used to mutually exclude threads' access to the VCF file
//  and overlapping genotypes for the current window.
pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;
// Used to mutually exclude appending Window_t* to the list of windows.
pthread_mutex_t windowListLock = PTHREAD_MUTEX_INITIALIZER;
// Used to mutually exclude access to the genome-wide IBS counts.
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

// All of the attributes used to window the genome.
typedef struct {
    VCFLocusParser_t* parser;
    HaplotypeEncoder_t* encoder;
    // Used to perform MDS on an asd matrix.
    RealSymEigen_t* eigen;
    // Dimension to project into.
    int k;

    // Our list of window pointers.
    klist_t(WindowPtr)* winList;

    int HAP_SIZE;
    int STEP_SIZE;
    int WINDOW_SIZE;
    int numSamples;

    // IMPORTANT!
    //  We create a queue of haplotypes to window the genome.
    //  winGeno is a WINDOW_SIZE-by-numSamples matrix. winStartIndex
    //  is the index of the first haplotype in winGeno for the current
    //  window. winEndIndex is the index of the last haplotype in winGeno
    //  for the current window. When we advance the window, we dequeue
    //  STEP_SIZE haplotypes. Then, as we form the current window, we 
    //  enqueue each new haplotype. We track the start loci of each 
    //  overlapping window in the array winStartLoci. We track the number
    //  of haplotypes already processed in overlapping windows with
    //  numHapsInOverlap.
    Genotype_t** winGeno;
    int winStartIndex;
    int winEndIndex;
    int* winStartLoci;
    int numHapsInOverlap;

    // The maximum number of basepairs a window can spread.
    int MAX_GAP;

    // Genome-wide IBS counts.
    IBS_t* globalAlleleCounts;
    // Number of haplotypes genome-wide.
    int globalNumHaps;
    // Number of loci genome-wide.
    int globalNumLoci;

    // The window number of the current window.
    int curWinNum;
    // The window number on the current chromosome.
    int curWinNumOnChrom;
    // The chromosome of the current window.
    kstring_t* curChrom;
} WindowRecord_t;

// Perform MDS on a window given its asd matrix.
// Accepts:
//  Window_t* window -> The current window.
//  RealSymEigen_t* eigen -> Used for MDS calculation.
//  double* asd -> The packed asd matrix.
//  int k -> The dimension to project down into.
// Returns: void.
void perform_mds_on_window(Window_t* window, RealSymEigen_t* eigen, double* asd, int k) {
    // If we drop a window due to it spanning to long of a region, then we do not perform MDS.
    if (window -> dropWindow) {
        LOG_WARNING("Drop window %d on %s because it exceeds gap threshold.\n", window -> winNumOnChrom, ks_str(window -> chromosome));
        window -> X = NULL;
        return;
    }
    // Create numSamples-by-k matrix.
    double** X = create_matrix(double, eigen -> N, k);
    // Compute MDS.
    int INFO = compute_classical_mds(eigen, asd, k, X);
    // If there was an error while computing MDS, destroy X and return.
    //  window -> X is NULL.
    if (INFO != 0) {
        LOG_WARNING("ASD matrix of window %d on %s had a NaN or rank less than %d. Could not calculate MDS.\n", window -> winNumOnChrom, ks_str(window -> chromosome), eigen -> N);
        destroy_matrix(double, X, eigen -> N);
        return;
    } else {
        LOG_INFO("Finished MDS for window %d on %s from %d to %d.\n", window -> winNum, ks_str(window -> chromosome), window -> startCoord, window -> endCoord);
    }
    window -> sigma = normalize_matrix(X, eigen -> N, k);
    window -> X = X;
}

// Reads in haplotypes and prepares the next window.
// Accepts:
//  WindowRecord_t* record -> The current state of the sliding window.
//  Genotype_t** threadGeno -> If single thread, pass NULL. When multiple 
//                              threads are used, each thread has their own
//                              matrix to hold the entire window for calculations.
// Returns: Window_t*, The current window.
Window_t* get_next_window(WindowRecord_t* record, Genotype_t** threadGeno) {

    // Create the new window and set bookkeeping values.
    Window_t* window = init_window();
    window -> winNum = record -> curWinNum++;
    window -> winNumOnChrom = record -> curWinNumOnChrom++;
    window -> chromosome = calloc(1, sizeof(kstring_t));
    kputs(ks_str(record -> curChrom), window -> chromosome);
    window -> numLoci = record -> numHapsInOverlap * record -> HAP_SIZE;

    LOG_INFO("Started ASD calculations for window %d.\n", window -> winNum);

    // Assume there will be another window on the same chromosome.
    bool isSameChrom = true;

    // Used for swapping.
    Genotype_t* temp;

    // If the window is the first one on the chromosome, we reset the queue.
    if (window -> winNumOnChrom == 1) {
        record -> winStartIndex = 0;
        record -> winEndIndex = 0;
    // Otherwise, we remove STEP_SIZE haplotypes from the queue.
    } else {
        record -> winStartIndex = (record -> winStartIndex + record -> STEP_SIZE) % record -> WINDOW_SIZE;
    }

    // While we have not reached that maximum number of haplotypes in the window, we
    //  have not reached EOF, and we did not reach the end of the chromosome, we
    //  read in the haplotypes that make up the current window.
    while(record -> numHapsInOverlap < record -> WINDOW_SIZE && isSameChrom) {
        // Get the next haplotype.
        isSameChrom = get_next_haplotype(record -> parser, record -> encoder, record -> HAP_SIZE);
        // If we are using a single thread, add haplotype to the end of the queue.
        if (threadGeno == NULL) {
            SWAP(record -> winGeno[record -> winEndIndex], record -> encoder -> genotypes, temp);
        // If we are using multiple threads, we need to add the haplotype to the end of the queue
        //  for future threads, but also, we need to copy its contents over to the operating thread's 
        //  matrix.
        } else {
            SWAP(record -> winGeno[record -> winEndIndex], record -> encoder -> genotypes, temp);
            for (int i = 0; i < record -> numSamples; i++)
                threadGeno[record -> numHapsInOverlap][i] = record -> winGeno[record -> winEndIndex][i];
        }
        record -> winEndIndex = (record -> winEndIndex + 1) % record -> WINDOW_SIZE;
        // If the haplotype is the start of an overlapping window, we save the start coordinate.
        if (record -> numHapsInOverlap % record -> STEP_SIZE == 0)
            record -> winStartLoci[(window -> winNumOnChrom + record -> numHapsInOverlap / record -> STEP_SIZE) % (record -> WINDOW_SIZE / record -> STEP_SIZE + 1)] = record -> encoder -> startCoord;
        window -> numLoci += record -> encoder -> numLoci;
        record -> globalNumLoci += record -> encoder -> numLoci;
        // We processed another haplotype.
        record -> numHapsInOverlap++;
    }

    // Set the start coordinate of the current window.
    window -> startCoord = record -> winStartLoci[window -> winNumOnChrom % (record -> WINDOW_SIZE / record -> STEP_SIZE + 1)];
    // Set the end coordinate of the current window.
    window -> endCoord = record -> encoder -> endCoord;
    // Drop window if it exceeds MAX_GAP.
    if (window -> endCoord - window -> startCoord > record -> MAX_GAP) {
        window -> dropWindow = true;
    }
    window -> numHaps = record -> numHapsInOverlap;
    if (isSameChrom) {
        record -> globalNumHaps += record -> STEP_SIZE;
        record -> numHapsInOverlap = record -> WINDOW_SIZE - record -> STEP_SIZE;
    // If we reached the end of the chromosome or file, reset chromosome and counters.
    } else {
        record -> globalNumHaps += window -> numLoci / record -> HAP_SIZE;
        ks_overwrite(ks_str(record -> parser -> nextChrom), record -> curChrom);
        record -> numHapsInOverlap = 0;
        record -> curWinNumOnChrom = 1;
    }

    return window;
}

// When the user specifies multiple threads, each
//  thread executes this method.
// Accepts:
//  void* arg -> Pointer to the created WindowRecord_t*.
// Returns: void.
void* sliding_window_multi_thread(void* arg) {

    WindowRecord_t* record = (WindowRecord_t*) arg;
    int numSamples = record -> numSamples;
    int HAP_SIZE = record -> HAP_SIZE;
    int STEP_SIZE = record -> STEP_SIZE;
    int k = record -> k;

    // Pointer to the current window.
    Window_t* window;

    // Create the genotype matrix to contain the entire current window.
    Genotype_t** threadGeno = create_matrix(geno, record -> WINDOW_SIZE, record -> numSamples);

    // IBS_t packed matrix to hold current window's IBS counts.
    IBS_t* winAlleleCounts = (IBS_t*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS_t));
    // IBS_t packed matrix to hold genome-wide related counts processed by thread.
    IBS_t* globalAlleleCounts = (IBS_t*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS_t));
    // Packed matrix to hold ASD values for MDS calculation.
    double* asd = (double*) calloc(PACKED_SIZE(record -> numSamples), sizeof(double));
    // Used to perform MDS calculations.
    RealSymEigen_t* eigen = init_real_sym_eigen(record -> numSamples);

    int start, numHapsInWin;
    bool isLastWinOnChrom;

    // While there are windows to process ...
    while (true) {

        // Get access to genome.
        pthread_mutex_lock(&genomeLock);
        // Once we have read the whole VCF, exit loop.
        if (record -> parser -> isEOF) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        // If the current window is not the first window on the chromsome, 
        //  we need to copy over the overlapping haplotypes from the previous window
        //  before we read in any additional haplotypes.
        if (record -> curWinNumOnChrom != 1) {
            start = (record -> winStartIndex + record -> STEP_SIZE) % record -> WINDOW_SIZE;
            for (int i = 0; i < record -> WINDOW_SIZE - record -> STEP_SIZE; i++) {
                for (int j = 0; j < record -> numSamples; j++)
                    threadGeno[i][j] = record -> winGeno[start][j];
                start = (start + 1) % record -> WINDOW_SIZE;
            }
        }
        // Get the rest of the window.
        window = get_next_window(record, threadGeno);
        // Set flag if the current window is the last window on the chromosome.
        isLastWinOnChrom = record -> curWinNumOnChrom == 1;
        pthread_mutex_unlock(&genomeLock);

        // Calculate the number of haplotypes in the window that we must process.
        numHapsInWin = (int) ceil((double) window -> numLoci / HAP_SIZE);

        // Compute ASD matrix for current window.
        process_window_multi_thread(threadGeno, winAlleleCounts, globalAlleleCounts, asd, window, numHapsInWin, isLastWinOnChrom, numSamples, STEP_SIZE);
        
        // Project samples into dimension k.
        perform_mds_on_window(window, eigen, asd, k);

        // Add the current window to the list.
        pthread_mutex_lock(&windowListLock);
        *kl_pushp(WindowPtr, record -> winList) = window;
        pthread_mutex_unlock(&windowListLock);

    }

    // Gain access to genome-wide IBS counts.
    pthread_mutex_lock(&globalLock);
    // Accumulate threads' genome-wide IBS contributions.
    for (int i = 0; i < record -> numSamples; i++)
        for (int j = i + 1; j < record -> numSamples; j++) 
            add_ibs(&(record -> globalAlleleCounts[PACKED_INDEX(i, j)]), &globalAlleleCounts[PACKED_INDEX(i, j)]);
    pthread_mutex_unlock(&globalLock);

    // Free all memory belonging to the thread.
    destroy_matrix(geno, threadGeno, record -> WINDOW_SIZE);
    free(winAlleleCounts);
    free(globalAlleleCounts);
    free(asd);
    destroy_real_sym_eigen(eigen);

    return NULL;

}

void sliding_window_single_thread(WindowRecord_t* record) {

    // Very similar to the multi-threaded version except, we must keep track 
    //  of the IBS counts in the first STEP_SIZE haplotypes. When we advance
    //  the window, we subtract stepAlleleCounts from winAlleleCounts to get
    //  the IBS counts of the haplotypes in the window overlap.

    Window_t* window;
    IBS_t* winAlleleCounts = (IBS_t*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS_t));
    IBS_t* stepAlleleCounts = (IBS_t*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS_t));
    double* asd = (double*) calloc(PACKED_SIZE(record -> numSamples), sizeof(double));

    int numHapsInWin;

    while (!(record -> parser -> isEOF)) {
        window = get_next_window(record, NULL);
        numHapsInWin = (int) ceil((double) window -> numLoci / record -> HAP_SIZE);
        process_window_single_thread(record -> winGeno, record -> winStartIndex, winAlleleCounts, stepAlleleCounts, record -> globalAlleleCounts, asd, window, numHapsInWin, window -> winNumOnChrom == 1, record -> curWinNumOnChrom == 1, record -> numSamples, record -> STEP_SIZE, record -> WINDOW_SIZE);
        perform_mds_on_window(window, record -> eigen, asd, record -> k);
        *kl_pushp(WindowPtr, record -> winList) = window;
    }

    free(winAlleleCounts);
    free(stepAlleleCounts);
    free(asd);
}

// Used to sort the list of windows when multiple threads are used.
// Accepts:
//  Window_t** windows -> The array of pointers to windows. Upon exit,
//                          the array will be sorted according to winNum.
//  int numWindows -> The number of windows in the array.
// Returns: void.
void sort_windows(Window_t** windows, int numWindows) {
    Window_t* temp;
    int j;
    // Using insertion sort since the windows should not be greatly out of order.
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

Window_t** sliding_window(VCFLocusParser_t* parser, HaplotypeEncoder_t* encoder, int k, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int MAX_GAP, int* numWindows) {
    
    // Set all attributes for the window record.
    WindowRecord_t* record = (WindowRecord_t*) calloc(1, sizeof(WindowRecord_t));
    record -> parser = parser;
    record -> encoder = encoder;
    record -> eigen = init_real_sym_eigen(encoder -> numSamples);
    record -> winList = kl_init(WindowPtr);

    record -> k = k;

    record -> HAP_SIZE = HAP_SIZE;
    record -> STEP_SIZE = STEP_SIZE;
    record -> WINDOW_SIZE = WINDOW_SIZE;
    record -> numSamples = encoder -> numSamples;
    record -> MAX_GAP = MAX_GAP;

    record -> winGeno = create_matrix(geno, WINDOW_SIZE, encoder -> numSamples);
    record -> numHapsInOverlap = 0;
    record -> winStartIndex = 0;
    record -> winEndIndex = 0;
    record -> winStartLoci = (int*) calloc(WINDOW_SIZE / STEP_SIZE + 1, sizeof(int));

    record -> globalAlleleCounts = (IBS_t*) calloc(PACKED_SIZE(record -> numSamples), sizeof(IBS_t));
    record -> globalNumHaps = 0;
    record -> globalNumLoci = 0;

    record -> curWinNum = 1;
    record -> curWinNumOnChrom = 1;
    record -> curChrom = calloc(1, sizeof(kstring_t));
    kputs(ks_str(parser -> nextChrom), record -> curChrom);

    // Our array of windows along the genome.
    Window_t** windows;

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
    
    // Create the genome-wide window.
    //  Will be the first window in the array.
    Window_t* global = init_window();
    global -> chromosome = calloc(1, sizeof(kstring_t));
    kputs("Global", global -> chromosome);
    global -> winNum = 0;
    global -> winNumOnChrom = 0;
    global -> startCoord = 0;
    global -> endCoord = 0;
    global -> sigma = -1;
    global -> numLoci = record -> globalNumLoci;
    global -> numHaps = record -> globalNumHaps;
    global -> saveIBS = true;
    global -> ibs = record -> globalAlleleCounts;

    double* globalASD = (double*) malloc(PACKED_SIZE(record -> numSamples) * sizeof(double));
    // Convert IBS counts to ASD.
    for (int i = 0; i < record -> numSamples; i++) 
        for (int j = i + 1; j < record -> numSamples; j++)
            globalASD[PACKED_INDEX(i, j)] = ibs_to_asd(record -> globalAlleleCounts[PACKED_INDEX(i, j)]);
    // Perfrom MDS.
    perform_mds_on_window(global, record -> eigen, globalASD, k);
    free(globalASD);

    // Convert the list of windows to an array.
    *numWindows = record -> curWinNum;
    windows = (Window_t**) calloc(*numWindows, sizeof(Window_t*));
    int i = 1;
    for (kliter_t(WindowPtr)* it = kl_begin(record -> winList); it != kl_end(record -> winList); it = kl_next(it)) {
        windows[i] = kl_val(it);
        i++;
    }
    windows[0] = global;
    // Sort the array of windows if multiple threads were used.
    if (NUM_THREADS > 1)
        sort_windows(windows, *numWindows);

    // Free all used memory.
    kl_destroy(WindowPtr, record -> winList);
    destroy_matrix(geno, record -> winGeno, WINDOW_SIZE);
    destroy_real_sym_eigen(record -> eigen);
    free(record -> winStartLoci);
    destroy_kstring(record -> curChrom);
    free(record);

    return windows;
}

// All attributes needed to just calculate the genome-wide window.
typedef struct {
    VCFLocusParser_t* parser;
    HaplotypeEncoder_t* encoder;
    IBS_t* alleleCounts;
    int numLoci;
    int numHaps;
    int HAP_SIZE;
} GlobalWindowRecord_t;

// When the user specifies multiple threads, each
//  thread executes this method.
// Accepts:
//  void* arg -> Pointer to the created GlobalWindowRecord_t*.
// Returns: void.
void* global_window_multi_thread(void* arg) {
    GlobalWindowRecord_t* record = (GlobalWindowRecord_t*) arg;
    int numSamples = record -> encoder -> numSamples;

    // IBS counts proccessed by a thread.
    IBS_t* alleleCounts = (IBS_t*) calloc(PACKED_SIZE(numSamples), sizeof(IBS_t));
    // We process each haplotype individually.
    Genotype_t* genotypes = (Genotype_t*) calloc(numSamples, sizeof(Genotype_t));
    // Used for swapping.
    Genotype_t* temp;

    while (true) {
        pthread_mutex_lock(&genomeLock);
        if (record -> parser -> isEOF) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        // Get the next haplotype.
        get_next_haplotype(record -> parser, record -> encoder, record -> HAP_SIZE);
        record -> numLoci += record -> encoder -> numLoci;
        record -> numHaps++;
        // Swap out haplotype so thread can calulate IBS.
        SWAP(genotypes, record -> encoder -> genotypes, temp);
        pthread_mutex_unlock(&genomeLock);
        // Process IBS of haplotype.
        pairwise_ibs(genotypes, alleleCounts, numSamples);
    }

    // Add counts to genome-wide IBS counts.
    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < numSamples; i++)
        for (int j = i + 1; j < numSamples; j++)
            add_ibs(&(record -> alleleCounts[PACKED_INDEX(i, j)]), &alleleCounts[PACKED_INDEX(i, j)]);
    pthread_mutex_unlock(&globalLock);

    free(alleleCounts);
    free(genotypes);
    return NULL;
}

Window_t* global_window(VCFLocusParser_t* parser, HaplotypeEncoder_t* encoder, int k, int HAP_SIZE, int NUM_THREADS) {

    // Genome-wide IBS counts.
    IBS_t* alleleCounts = (IBS_t*) calloc(PACKED_SIZE(encoder -> numSamples), sizeof(IBS_t));
    // Used to hold pairwise ASD values for MDS.
    double* asd = (double*) calloc(PACKED_SIZE(encoder -> numSamples), sizeof(double));

    // Create global window.
    Window_t* window = init_window();
    window -> winNum = 0;
    window -> winNumOnChrom = 0;
    window -> chromosome = calloc(1, sizeof(kstring_t));
    kputs("Global", window -> chromosome);
    window -> startCoord = 0;
    window -> endCoord = 0;

    // If we are using one thread, directly calculate genome-wide.
    if (NUM_THREADS == 1) {
        while (!parser -> isEOF) {
            get_next_haplotype(parser, encoder, HAP_SIZE);
            window -> numLoci += encoder -> numLoci;
            window -> numHaps++;
            pairwise_ibs(encoder -> genotypes, alleleCounts, encoder -> numSamples);
        }
    } else {
        GlobalWindowRecord_t* record = (GlobalWindowRecord_t*) calloc(1, sizeof(GlobalWindowRecord_t));
        record -> parser = parser;
        record -> encoder = encoder;
        record -> alleleCounts = alleleCounts;
        record -> numLoci = 0;
        record -> HAP_SIZE = HAP_SIZE;
        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        // Create threads.
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_create(&threads[i], NULL, global_window_multi_thread, (void*) record);
        global_window_multi_thread((void*) record);
        // Destroy threads.
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);
        window -> numLoci = record -> numLoci;
        window -> numHaps = record -> numHaps;
        free(threads);
        free(record);
    }

    LOG_INFO("Finished genome-wide ASD calculations. Performing MDS ...\n");

    // Convert counts to ASD.
    for (int i = 0; i < encoder -> numSamples; i++)
        for (int j = i + 1; j < encoder -> numSamples; j++)
            asd[PACKED_INDEX(i, j)] = ibs_to_asd(alleleCounts[PACKED_INDEX(i, j)]);

    // Perfrom MDS on global.
    RealSymEigen_t* eigen = init_real_sym_eigen(encoder -> numSamples);
    perform_mds_on_window(window, eigen, asd, k);
    destroy_real_sym_eigen(eigen);

    // Save ibs matrix for global.
    window -> saveIBS = true;
    window -> ibs = alleleCounts;

    free(asd);

    return window;
}

/*
void print_window_info(Window_t* window, int n, int k) {
    printf("Window Number: %d\n", window -> winNum);
    printf("Chromosome: %s\n", ks_str(window -> chromosome));
    printf("Window Number on Chromosome: %d\n", window -> winNumOnChrom);
    printf("Start Position: %d\n", window -> startCoord);
    printf("End Position: %d\n", window -> endCoord);
    printf("Number of Loci: %d\n", window -> numLoci);
    printf("X = \n");
    for (int i = 0; i < n && window -> X != NULL; i++) {
        for (int j = 0; j < k; j++)
            printf("%5lf\t", window -> X[i][j]);
        printf("\n");
    }
    printf("\n\n");
}

int main() {

    int k = 2;
    int HAP_SIZE = 1;
    int STEP_SIZE = 2; 
    int WINDOW_SIZE = 3;
    int NUM_THREADS = 2;
    int numWindows;

    VCFLocusParser_t* parser = init_vcf_locus_parser("./data/sliding_window_test2.vcf.gz", NULL, false, 0, 1, true);
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(parser -> numSamples);

    Window_t** windows = sliding_window(parser, encoder, k, HAP_SIZE, STEP_SIZE, WINDOW_SIZE, NUM_THREADS, &numWindows);

    for (int i = 0; i < numWindows; i++) {
        print_window_info(windows[i], encoder -> numSamples, k);
        destroy_window(windows[i], encoder -> numSamples);
    }
    free(windows);

    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);

}
*/