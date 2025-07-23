// File: BlockList.h
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#ifndef _BLOCK_LIST_H_
#define _BLOCK_LIST_H_

#include <stdbool.h>

// Our structure that tracks the counts of loci
//  with 0, 1, and 2 alleles IBS between samples.
typedef struct {
    unsigned int ibs0;
    unsigned int ibs1;
    unsigned int ibs2;
} IBS_t;

// A node in the BlockList.
typedef struct Block {
    // Allele counts and lower dimensional representation for the block.
    IBS_t* alleleCounts;
    double** X;

    // Procrustes t-statistic, t-statistic if block is excluded, and effective rank of points.
    double procrustesT;
    double excludedT;
    double effectRank;

    // Block attributes.
    int blockNum;
    int blockNumOnChrom;
    char* chrom;
    int startCoordinate;
    int endCoordinate;
    int numHaps;
    bool isDropped;
    int numSamples;
    // The next block in the list.
    struct Block* next;
} Block_t;

// A list of blocks.
typedef struct BlockList {
    // Global counts.
    IBS_t* alleleCounts;
    double** X;
    double effectRank;

    // If jackknife is computed
    double procrustesT;
    double stdDev;
    double pvalue;

    // Global attributes.
    int numSamples;
    int numBlocks;
    int numHaps;
    Block_t* head;
    Block_t* tail;
} BlockList_t;

// Creates a list of blocks.
// Accepts:
//  int numSamples -> The number of samples.
// Returns: An empty block list.
BlockList_t* init_block_list(int numSamples);

// Creates a block.
// Acccepts:
//  char* chrom -> The chromosome the block is on.
//  int startCoordinate -> The coordinate of the first record in the block.
//  int numSamples -> The number of samples.
// Returns: Block_t*, the new block.
Block_t* init_block(char* chrom, int startCoordinate, int numSamples);

// Adds a block to the block list.
// Accepts:
//  BlockList_t* blockList -> The list to add the block to.
//  Block_t* block -> Adds the block to the end of the list.
// Returns: void.
void append_block(BlockList_t* blockList, Block_t* block);

// Frees the memory occupied by the block list.
// Accepts:
//  BlockList_t* blockList -> The block list to destroy.
// Returns: void.
void destroy_block_list(BlockList_t* blockList);

#endif