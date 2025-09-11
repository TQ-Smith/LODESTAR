// File: LODESTAR.c
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main LODESTAR operations.

#include "LODESTAR.h"
#include <pthread.h>
#include <stdio.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

// Used to mutually exclude threads' access to the VCF file
//  and overlapping genotypes for the current window.
pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;

// Used to mutually exclude access to the global list.
pthread_mutex_t listLock = PTHREAD_MUTEX_INITIALIZER;

// Used to mutually exclude access to the global counts.
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

// Saves a few lines to swap values.
#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

// A list of haplotypes within a block.
typedef struct BlockGenotypes {
    Genotype_t* genotypes;
    struct BlockGenotypes* next;
} BlockGenotypes_t;

// What is given to each thread to block the genome.
typedef struct BlockCounts {

    // A simple linked list to add haplotypes to the block.
    //  We save a little time by not reallocating memory for each block.
    BlockGenotypes_t* head;
    BlockGenotypes_t* tail;
    int numHapsInGenotypeList;
    int numHaps;

    // Protected by mutex for reading in.
    VCFLocusParser_t* vcfFile;
    HaplotypeEncoder_t* encoder;
    int numSamples;
    int blockSize;
    int haplotypeSize;
    int dropThreshold;

    // Protected by mutex for adding values to global.
    BlockList_t* globalList;
    int* blockNum;

} BlockCounts_t;

// A simple method to add a new haplotype to the linked list within a block.
void add_haplotype_to_block(BlockCounts_t* blockCounts) {
    BlockGenotypes_t* temp = calloc(1, sizeof(BlockGenotypes_t));
    temp -> genotypes = calloc(blockCounts -> numSamples, sizeof(Genotype_t));
    temp -> next = NULL;
    if (blockCounts -> numHapsInGenotypeList == 0) {
        blockCounts -> head = temp;
        blockCounts -> tail = temp;
    } else {
        blockCounts -> tail -> next = temp;
        blockCounts -> tail = temp;
    }
    blockCounts -> numHapsInGenotypeList++;
}

// Frees the linked list of haplotypes within a block.
void destroy_block_counts(BlockCounts_t* blockCounts) {
    if (blockCounts -> head != NULL) {
        BlockGenotypes_t* temp = NULL;
        for (int i = 0; i < blockCounts -> numHapsInGenotypeList; i++) {
            temp = blockCounts -> head;
            blockCounts -> head = blockCounts -> head -> next;
            free(temp -> genotypes);
            free(temp);
        }
    }
    if (blockCounts != NULL)
        free(blockCounts);
}

// For each thread, we read in a block of haplotypes and compute ASD.
void* partition(void* arg) {
    BlockCounts_t* blockCounts = (BlockCounts_t*) arg;

    Genotype_t* genoTemp;
    BlockGenotypes_t* blockTemp;

    // Points to the new block we add to the global list.
    Block_t* block;

    // Holds global counts.
    IBS_t* globalCounts = calloc(PACKED_SIZE(blockCounts -> numSamples), sizeof(IBS_t));

    while (true) {

        pthread_mutex_lock(&genomeLock);
        if (isEOF(blockCounts -> vcfFile)) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        
        // Reset everything for next block.
        block = init_block(blockCounts -> vcfFile -> nextChrom, blockCounts -> vcfFile -> nextCoord, blockCounts -> numSamples);
        blockCounts -> numHaps = 0;
        blockTemp = blockCounts -> head;
        bool isOnSameChrom = true;
        int endOfBlock = ((int) ((blockCounts -> vcfFile -> nextCoord - 1) / (double) blockCounts -> blockSize) + 1) * blockCounts -> blockSize;
        
        // Read in the next block.
        while (isOnSameChrom && blockCounts -> vcfFile -> nextCoord <= endOfBlock) {
            isOnSameChrom = get_next_haplotype(blockCounts -> vcfFile, blockCounts -> encoder, blockCounts -> haplotypeSize);
            // If we need more haplotypes, then we allocate the necessary memory.
            if (blockCounts -> numHaps == blockCounts -> numHapsInGenotypeList) {
                add_haplotype_to_block(blockCounts);
                blockTemp = blockCounts -> tail;
            }
            // Insert the newly read in haplotype.
            SWAP(blockTemp -> genotypes, blockCounts -> encoder -> genotypes, genoTemp);
            blockTemp = blockTemp -> next;
            blockCounts -> numHaps++;
            block -> numHaps++;
            block -> numLoci += blockCounts -> encoder -> numLoci;
        }
        // Finish the block.
        block -> endCoordinate = blockCounts -> encoder -> endCoord;
        block -> blockNum = *(blockCounts -> blockNum);
        *(blockCounts -> blockNum) = *(blockCounts -> blockNum) + 1;
        if (block -> numHaps < blockCounts -> dropThreshold)
            block -> isDropped = true;
        pthread_mutex_unlock(&genomeLock);

        // Calculate IBS for all haplotypes within the block and accumulate global if block is not dropped.
        if (!block -> isDropped) {
            block -> alleleCounts = calloc(PACKED_SIZE(blockCounts -> numSamples), sizeof(IBS_t));
            BlockGenotypes_t* locus = blockCounts -> head;
            for (int l = 0; l < blockCounts -> numHaps; l++) {
                for (int i = 0; i < blockCounts -> numSamples; i++) {
                    for (int j = i + 1; j < blockCounts -> numSamples; j++) {
                        if (locus -> genotypes[i].left != MISSING || locus -> genotypes[j].right != MISSING) {
                            int numSharedAlleles = num_shared_alleles(locus -> genotypes[i], locus -> genotypes[j]);
                            increment_ibs_value(&(block -> alleleCounts[PACKED_INDEX(i, j)]), numSharedAlleles);
                            increment_ibs_value(&(globalCounts[PACKED_INDEX(i, j)]), numSharedAlleles);
                        }
                    }
                }
                locus = locus -> next;
            }
        }

        // Append the block to the global list.
        pthread_mutex_lock(&listLock);
        if (block -> isDropped)
            fprintf(stderr, "Block on %s from %d to %d contains %d haplotypes. Dropped block.\n", block -> chrom, block -> startCoordinate, block -> endCoordinate, block -> numHaps);
        else
            fprintf(stderr, "Finished IBS for block number %d on %s from %d to %d.\n", block -> blockNum, block -> chrom, block -> startCoordinate, block -> endCoordinate);
        append_block(blockCounts -> globalList, block);
        pthread_mutex_unlock(&listLock);
    }

    // Accumulate global counts.
    pthread_mutex_lock(&globalLock);
    for (int i = 0; i < blockCounts -> numSamples; i++)
        for (int j = i + 1; j < blockCounts -> numSamples; j++)
            add_ibs(&(blockCounts -> globalList -> alleleCounts[PACKED_INDEX(i, j)]), &(globalCounts[PACKED_INDEX(i, j)]));
    pthread_mutex_unlock(&globalLock);

    destroy_block_counts(blockCounts);
    free(globalCounts);

    return NULL;
}

BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int dropThreshold, int NUM_THREADS) {
    
    BlockList_t* globalList = init_block_list(numSamples);
    globalList -> alleleCounts = calloc(PACKED_SIZE(numSamples), sizeof(IBS_t));
    int* blockNum = calloc(1, sizeof(int));
    *blockNum = 1;
    
    if (NUM_THREADS == 1) {
        BlockCounts_t* blockCounts = calloc(1, sizeof(BlockCounts_t));
        blockCounts -> vcfFile = vcfFile;
        blockCounts -> encoder = encoder;
        blockCounts -> blockSize = blockSize;
        blockCounts -> numSamples = numSamples;
        blockCounts -> haplotypeSize = haplotypeSize;
        blockCounts -> dropThreshold = dropThreshold;
        blockCounts -> globalList = globalList;
        blockCounts -> blockNum = blockNum;
        partition((void*) blockCounts);
    } else {
        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++) {
            BlockCounts_t* blockCounts = calloc(1, sizeof(BlockCounts_t));
            blockCounts -> vcfFile = vcfFile;
            blockCounts -> encoder = encoder;
            blockCounts -> blockSize = blockSize;
            blockCounts -> numSamples = numSamples;
            blockCounts -> haplotypeSize = haplotypeSize;
            blockCounts -> dropThreshold = dropThreshold;
            blockCounts -> globalList = globalList;
            blockCounts -> blockNum = blockNum;
            pthread_create(&threads[i], NULL, partition, (void*) blockCounts);
        }
        BlockCounts_t* blockCounts = calloc(1, sizeof(BlockCounts_t));
        blockCounts -> vcfFile = vcfFile;
        blockCounts -> encoder = encoder;
        blockCounts -> blockSize = blockSize;
        blockCounts -> numSamples = numSamples;
        blockCounts -> haplotypeSize = haplotypeSize;
        blockCounts -> dropThreshold = dropThreshold;
        blockCounts -> globalList = globalList;
        blockCounts -> blockNum = blockNum;
        partition((void*) blockCounts);
        // Wait for all threads to finish.
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);

        // Selection sort on block number.
        //  Didn't want to think about this.
        //  https://www.geeksforgeeks.org/dsa/iterative-selection-sort-for-linked-list/
        Block_t* sortedHead = NULL;
        Block_t* sortedTail = NULL;
        while (globalList -> head != NULL) {
            Block_t* min = globalList -> head;
            Block_t* prevMin = NULL;
            Block_t* current = globalList -> head;
            Block_t* prev = NULL;
            while (current != NULL) {
                if (current -> blockNum < min -> blockNum) {
                    min = current;
                    prevMin = prev;
                }
                prev = current;
                current = current -> next;
            }
            if (min == globalList -> head) {
                globalList -> head = globalList -> head -> next;
            } else {
                prevMin -> next = min -> next;
            }
            if (sortedHead == NULL) {
                sortedHead = min;
                sortedTail = min;
            } else {
                sortedTail -> next = min;
                sortedTail = min;
            }
        }
        globalList -> head = sortedHead;
    }
    free(blockNum);

    // Assign blockNumOnChrom and count global number of haps.
    int blockNumOnChrom = 1;
    globalList -> numHaps = 0;
    for (Block_t* temp = globalList -> head; temp -> next != NULL; temp = temp -> next) {
        temp -> blockNumOnChrom = blockNumOnChrom;
        if (strcmp(temp -> chrom, temp -> next -> chrom) != 0)
            blockNumOnChrom = 1;
        else 
            blockNumOnChrom++;
        if (!temp -> isDropped) {
            globalList -> numHaps += temp -> numHaps;
            globalList -> numLoci += temp -> numLoci;
        }
    }
    globalList -> tail -> blockNumOnChrom = blockNumOnChrom;

    return globalList;
}

typedef struct BlockProcrustes {
    BlockList_t* globalList;
    double** y;
    double* y0;
    int k;
    // Points to the current block in globalList for the next thread to operate on.
    Block_t** current;
    // For the bootstrap.
    int numReps;
    int sampleSize;
    int* currentReplicate;
} BlockProcrustes_t;

// Randomly sample a non-dropped block.
Block_t* get_random_block(gsl_rng* r, BlockList_t* globalList) {
    Block_t* temp = globalList -> head;
    while (true) {
        int blockNum = globalList -> numBlocks * (double) gsl_rng_uniform(r);
        for (int i = 0; i < blockNum; i++)
            temp = temp -> next;
        if (!temp -> isDropped)
            break;
    }
    return temp;
}

void* procrustes_bootstrap(void* arg) {
    BlockProcrustes_t* blockProcrustes = (BlockProcrustes_t*) arg;

    // Allocate all required memory.
    double* asdBlock = calloc(PACKED_SIZE(blockProcrustes -> globalList -> numSamples), sizeof(double));
    IBS_t* ibsBlock = calloc(PACKED_SIZE(blockProcrustes -> globalList -> numSamples), sizeof(IBS_t));
    double** bootX = init_matrix(blockProcrustes -> globalList -> numSamples, blockProcrustes -> k);
    RealSymEigen_t* eigen = init_real_sym_eigen(blockProcrustes -> globalList -> numSamples);

    // The current block to operate on.
    Block_t* current = NULL;

    while (true) {

        // Reuse list lock to get the next block.
        pthread_mutex_lock(&listLock);
        if (*(blockProcrustes -> current) == NULL) {
            pthread_mutex_unlock(&listLock);
            break;
        }
        current = *(blockProcrustes -> current);
        *(blockProcrustes -> current) = (*(blockProcrustes -> current)) -> next;
        if (current -> isDropped)
            fprintf(stderr, "Block on %s from %d to %d was dropped. Skipping.\n", current -> chrom, current -> startCoordinate, current -> endCoordinate);
        else
            fprintf(stderr, "Performing MDS and Procrustes for block number %d on %s from %d to %d.\n", current -> blockNum, current -> chrom, current -> startCoordinate, current -> endCoordinate);
        // If we reached the end of the list.
        if (current -> next == NULL && blockProcrustes -> numReps == 0)
            fprintf(stderr, "\nProcrustes finished ...\n\n");
        else if (current -> next == NULL)
            fprintf(stderr, "\nProcrustes finished. Starting bootstrap ...\n\n");
        pthread_mutex_unlock(&listLock);

        if (!current -> isDropped) {
            // Convert to ASD first.
            for (int i = 0; i < blockProcrustes -> globalList -> numSamples; i++)
                for (int j = i + 1; j < blockProcrustes -> globalList -> numSamples; j++)
                    asdBlock[PACKED_INDEX(i, j)] = ibs_to_asd(current -> alleleCounts[PACKED_INDEX(i, j)]);

            // Perform MDS and Procrustes.
            double** X = init_matrix(eigen -> N, blockProcrustes -> k);
            current -> varCapt = compute_classical_mds(eigen, asdBlock, blockProcrustes -> k, X);
            // I am doing this to be safe. Results of cMDS are already centered.
            normalize_matrix(X, blockProcrustes -> globalList -> numSamples, blockProcrustes -> k);

            // In case MDS did not converge, treat the block as having no effect.
            if (current -> varCapt == -1) {
                current -> X = NULL;
                current -> procrustesT = 0;
                destroy_matrix(X, eigen -> N);
            // If we are comparing against the global.
            } else if (blockProcrustes -> y == NULL) {
                current -> X = X;
                current -> procrustesT = procrustes_statistic(current -> X, NULL, blockProcrustes -> globalList -> X, NULL, eigen, eigen -> N, blockProcrustes -> k, true);
            // If we are comparing against user defined coordinates.
            } else {
                current -> X = X;
                current -> procrustesT = procrustes_statistic(current -> X, NULL, blockProcrustes -> y, blockProcrustes -> y0, eigen, eigen -> N, blockProcrustes -> k, true);
            }
        }
    }

    // Create RNG.
    gsl_rng_env_setup();
    const gsl_rng_type* T = gsl_rng_default;
    gsl_rng* r = gsl_rng_alloc(T);
    // Seed each thread differently.
    gsl_rng_set(r, time(NULL) + (long) arg);

    // Start the bootstrap.
    // Only execute 
    while (blockProcrustes -> numReps != 0) {
        double bootStrappedT;

        // Create our random replicate.
        for (int s = 0; s < blockProcrustes -> sampleSize; s++) {
            Block_t* temp = get_random_block(r, blockProcrustes -> globalList);
            for (int i = 0; i < blockProcrustes -> globalList -> numSamples; i++) {
                for (int j = i + 1; j < blockProcrustes -> globalList -> numSamples; j++) {
                    if (s == 0)
                        ibsBlock[PACKED_INDEX(i, j)] = temp -> alleleCounts[PACKED_INDEX(i, j)];
                    else
                        add_ibs(&(ibsBlock[PACKED_INDEX(i, j)]), &(temp -> alleleCounts[PACKED_INDEX(i, j)]));
                    if (s == blockProcrustes -> sampleSize - 1)
                        asdBlock[PACKED_INDEX(i, j)] = ibs_to_asd(temp -> alleleCounts[PACKED_INDEX(i, j)]);
                }
            }
        }

        // Perfrom MDS.
        double effectiveRank = compute_classical_mds(eigen, asdBlock, blockProcrustes -> k, bootX);

        // If MDS does not converge, then we do not count the sample.
        if (effectiveRank == -1)
            continue;

        normalize_matrix(bootX, blockProcrustes -> globalList -> numSamples, blockProcrustes -> k);

        // Calculate our Procrustes statistic.
        if (blockProcrustes -> y == NULL) 
            bootStrappedT = procrustes_statistic(bootX, NULL, blockProcrustes -> globalList -> X, NULL, eigen, eigen -> N, blockProcrustes -> k, false);
        else 
            bootStrappedT = procrustes_statistic(bootX, NULL, blockProcrustes -> y, blockProcrustes -> y0, eigen, eigen -> N, blockProcrustes -> k, false);
        

        pthread_mutex_lock(&genomeLock);
        if (*(blockProcrustes -> currentReplicate) == blockProcrustes -> numReps) {
            pthread_mutex_unlock(&genomeLock);
            break;
        }
        blockProcrustes -> globalList -> samplingDistribution[*(blockProcrustes -> currentReplicate)] = bootStrappedT;
        *(blockProcrustes -> currentReplicate) += 1;
        fprintf(stderr, "Completed Replicate %d of the bootstrap.\n", *(blockProcrustes -> currentReplicate));
        pthread_mutex_unlock(&genomeLock);
    }

    destroy_real_sym_eigen(eigen);
    destroy_matrix(bootX, eigen -> N);
    free(asdBlock);
    free(ibsBlock);
    free(blockProcrustes);
    gsl_rng_free(r);

    return NULL;
}

void procrustes(BlockList_t* globalList, double** y, double* y0, int k, int NUM_THREADS, int numReps, int sampleSize) {
    
    int* currentReplicate = calloc(1, sizeof(int));
    *currentReplicate = 0;
    globalList -> samplingDistribution = calloc(numReps, sizeof(double));
    Block_t* current = calloc(1, sizeof(Block_t*));
    current = globalList -> head;

    // Compute Procrustes statistic for each block and bootstrap.
    if (NUM_THREADS == 1) {
        BlockProcrustes_t* blockProcrustes = calloc(1, sizeof(BlockProcrustes_t));
        blockProcrustes -> globalList = globalList;
        blockProcrustes -> y = y;
        blockProcrustes -> y0 = y0;
        blockProcrustes -> k = k;
        blockProcrustes -> current = &current;
        blockProcrustes -> numReps = numReps;
        blockProcrustes -> sampleSize = sampleSize;
        blockProcrustes -> currentReplicate = currentReplicate;
        procrustes_bootstrap((void*) blockProcrustes);
    } else {
        pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
        for (int i = 0; i < NUM_THREADS - 1; i++) {
            BlockProcrustes_t* blockProcrustes = calloc(1, sizeof(BlockProcrustes_t));
            blockProcrustes -> globalList = globalList;
            blockProcrustes -> y = y;
            blockProcrustes -> y0 = y0;
            blockProcrustes -> k = k;
            blockProcrustes -> current = &current;
            blockProcrustes -> numReps = numReps;
            blockProcrustes -> sampleSize = sampleSize;
            blockProcrustes -> currentReplicate = currentReplicate;
            pthread_create(&threads[i], NULL, procrustes_bootstrap, (void*) blockProcrustes);
        }
        BlockProcrustes_t* blockProcrustes = calloc(1, sizeof(BlockProcrustes_t));
        blockProcrustes -> globalList = globalList;
        blockProcrustes -> y = y;
        blockProcrustes -> y0 = y0;
        blockProcrustes -> k = k;
        blockProcrustes -> current = &current;
        blockProcrustes -> numReps = numReps;
        blockProcrustes -> sampleSize = sampleSize;
        blockProcrustes -> currentReplicate = currentReplicate;
        procrustes_bootstrap((void*) blockProcrustes);
        // Wait for all threads to finish.
        for (int i = 0; i < NUM_THREADS - 1; i++)
            pthread_join(threads[i], NULL);
        free(threads);
    }
    free(currentReplicate);
    free(current);

    // Compute p-values. We are doing this the lazy way.
    for (Block_t* temp = globalList -> head; temp != NULL; temp = temp -> next) {
        int numGreater = 1;
        for (int j = 0; j < numReps; j++) {
            if (temp -> procrustesT <= globalList -> samplingDistribution[j])
                numGreater++;
        }
        temp -> pvalue = numGreater / (double) (numReps + 1);
    }
    // Global p-value.
    if (y != NULL) {
        int numGreater = 1;
        for (int j = 0; j < numReps; j++) {
            if (globalList -> procrustesT <= globalList -> samplingDistribution[j])
                numGreater++;
        }
        globalList -> pvalue = numGreater / (double) (numReps + 1);
    }
}