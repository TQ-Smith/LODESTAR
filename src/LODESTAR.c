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

// Used to mutually exclude access to the genome-wide IBS counts.
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


        pthread_mutex_lock(&globalLock);
        append_block(blockCounts -> globalList, block);
        pthread_mutex_unlock(&globalLock);
    }
    destroy_block_counts(blockCounts);
    return NULL;
}

BlockList_t* block_allele_sharing(VCFLocusParser_t* vcfFile, HaplotypeEncoder_t* encoder, int numSamples, int blockSize, int haplotypeSize, int NUM_THREADS) {
    
    BlockList_t* globalList = init_block_list(numSamples);
    BlockList_t* sortedGlobalList = calloc(1, sizeof(BlockList_t));
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

        sortedGlobalList -> numSamples = globalList -> numSamples;
        sortedGlobalList -> numHaps = globalList -> numHaps;
        sortedGlobalList -> alleleCounts = globalList -> alleleCounts;
        sortedGlobalList -> X = globalList -> X;
        Block_t* temp1, temp2;
        // Selection sort linked list of blocks.
        for (int i = 1; i <= globalList -> numBlocks - 1; i++) {
            for (Block_t* j = globalList; j -> next -> blockNum != i; j = j -> next) {
                temp1 = j;
                temp2 = j -> next;
            }
            temp1 -> next = temp1 -> next -> next;
            if (i == 1) {
                sortedGlobalList -> head = temp2;
                sortedGlobalList -> tail = temp2;
            } else if (i == globalList -> numBlocks - 1) {
                sortedGlobalList -> tail = globalList -> head;
            } else {
                sortedGlobalList -> tail -> next = temp2;
                sortedGlobalList -> tail = temp2;
            }
        }
    }
    free(globalList);
    free(blockNum);

    int blockNumOnChrom = 1;
    sortedGlobalList -> head -> numBlockOnChrom = 1;
    for (temp1 = sortedGlobalList -> head -> next; temp1 != NULL; temp1 = temp1 -> next) {
        
    }

    return sortedGlobalList;
}