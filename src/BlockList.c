// File: BlockList.c
// Date: 7 July 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define the list of blocks along the genome.

#include "BlockList.h"
#include <string.h>
#include <stdlib.h>

Block_t* init_block(char* chrom, int startCoordinate, int numSamples) {
    Block_t* block = calloc(1, sizeof(Block_t));
    block -> chrom = strdup(chrom);
    block -> startCoordinate = startCoordinate;
    block -> numHaps = 0;
    block -> numSamples = numSamples;
    block -> alleleCounts = NULL;
    block -> X = NULL;
    block -> next = NULL;
    return block;
}

BlockList_t* init_block_list(int numSamples) {
    BlockList_t* blockList = calloc(1, sizeof(BlockList_t));
    blockList -> numSamples = numSamples;
    blockList -> head = NULL;
    blockList -> tail = NULL;
    blockList -> numHaps = 0;
    blockList -> alleleCounts = NULL;
    blockList -> X = NULL;
    return blockList;
}

void append_block(BlockList_t* blockList, Block_t* block) {
    if (blockList -> numBlocks == 0) {
        blockList -> head = block;
        blockList -> tail = block;
    } else {
        blockList -> tail -> next = block;
        blockList -> tail = block;
    }
    blockList -> numBlocks++;
}

void destroy_block(Block_t* block) {
    if (block -> X != NULL) {
        for (int i = 0; i < block -> numSamples; i++) {
            free(block -> X[i]);
        }
        free(block -> X);
    }
    if (block -> alleleCounts != NULL)
        free(block -> alleleCounts);
    free(block -> chrom);
    free(block);
}

void destroy_block_list(BlockList_t* blockList) {
    if (blockList -> X != NULL) {
        for (int i = 0; i < blockList -> numSamples; i++) {
            free(blockList -> X[i]);
        }
        free(blockList -> X);
    }
    if (blockList -> alleleCounts != NULL)
        free(blockList -> alleleCounts);
    Block_t* temp = NULL;
    for (int i = 0; i < blockList -> numBlocks; i++) {
        temp = blockList -> head;
        blockList -> head = blockList -> head -> next;
        destroy_block(temp);
    }
    free(blockList);
}