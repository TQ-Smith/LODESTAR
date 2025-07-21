// File: LODESTAR.c
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Main LODESTAR operations.

#include "LODESTAR.h"
#include <pthread.h>
#include <stdio.h>

// Used to mutually exclude threads' access to the VCF file
//  and overlapping genotypes for the current window.
pthread_mutex_t genomeLock = PTHREAD_MUTEX_INITIALIZER;

// Used to mutually exclude access to the global list.
pthread_mutex_t listLock = PTHREAD_MUTEX_INITIALIZER;

// Used to mutually exclude access to the global counts.
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

// The number of elements in the upper triangle of a symmetric matrix.
#define PACKED_SIZE(N) ((N * (N + 1)) / 2)

// We store the upper triangle of a symmetric matrix in a one-dimensional array,
//  This is known as packed storage. Element Aij -> a_{i + j * (j + 1) / 2}
#define PACKED_INDEX(i, j) (i + j * (j + 1) / 2)

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
        }
        // Finish the block.
        block -> endCoordinate = blockCounts -> encoder -> endCoord;
        block -> blockNum = *(blockCounts -> blockNum);
        *(blockCounts -> blockNum) = *(blockCounts -> blockNum) + 1;
        pthread_mutex_unlock(&genomeLock);

        // Calculate IBS for all haplotypes within the block and accumulate global.
        block -> alleleCounts = calloc(PACKED_SIZE(blockCounts -> numSamples), sizeof(IBS_t));
        BlockGenotypes_t* locus = blockCounts -> head;
        for (int l = 0; l < blockCounts -> numHaps; l++) {
            for (int i = 0; i < blockCounts -> numSamples; i++) {
                for (int j = i + 1; j < blockCounts -> numSamples; j++) {
                    if (locus -> genotypes[i].left != MISSING && locus -> genotypes[j].right != MISSING) {
                        int numSharedAlleles = num_shared_alleles(locus -> genotypes[i], locus -> genotypes[j]);
                        increment_ibs_value(&(block -> alleleCounts[PACKED_INDEX(i, j)]), numSharedAlleles);
                        increment_ibs_value(&(globalCounts[PACKED_INDEX(i, j)]), numSharedAlleles);
                    }
                }
            }
            locus = locus -> next;
        }

        // Append the block to the global list.
        pthread_mutex_lock(&listLock);
        fprintf(stderr, "Finished block number %d on %s from %d to %d.\n", block -> blockNum, block -> chrom, block -> startCoordinate, block -> endCoordinate);
        append_block(blockCounts -> globalList, block);
        pthread_mutex_unlock(&listLock);
    }

    // Accumulate global counts.
    pthread_mutex_lock(&listLock);
    for (int i = 0; i < blockCounts -> numSamples; i++)
        for (int j = i + 1; j < blockCounts -> numSamples; j++)
            add_ibs(&(blockCounts -> globalList -> alleleCounts[PACKED_INDEX(i, j)]), &(globalCounts[PACKED_INDEX(i, j)]));
    pthread_mutex_unlock(&listLock);

    destroy_block_counts(blockCounts);
    free(globalCounts);

    return NULL;
}

BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int NUM_THREADS) {
    
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
        blockCounts -> globalList = globalList;
        blockCounts -> blockNum = blockNum;
        partition((void*) blockCounts);
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

    // Assign blockNumOnChrom.
    int blockNumOnChrom = 1;
    for (Block_t* temp = globalList -> head; temp -> next != NULL; temp = temp -> next) {
        temp -> blockNumOnChrom = blockNumOnChrom;
        if (strcmp(temp -> chrom, temp -> next -> chrom) != 0)
            blockNumOnChrom = 1;
        else 
            blockNumOnChrom++;
    }
    globalList -> tail -> blockNumOnChrom = blockNumOnChrom;

    return globalList;
}

typedef struct BlockProcrustes {

} BlockProcrustes_t;

void procrustes(BlockList_t* globalList, double** y, double* y0, int dropThreshold, int NUM_THREADS) {

}