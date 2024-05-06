
// File: RegionSet.h
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Used to define a set of genomic coordinates.

#ifndef _REGION_SET_H_
#define _REGION_SET_H_

#include "../lib/kstring.h"
#include "../lib/khash.h"
#include "stdbool.h"

// A region on a chromosome.
//  Used as a node in a SORTED linked-list.
typedef struct Region {
    unsigned int startLocus;
    unsigned int endLocus;
    // Points to next interval on the same chromosome such that
    //  startLocus <= next -> startLocus.
    struct Region* next;
} Region_t;

// Initialize a string hashtable that will associate chromosome name
//  with a sorted list of intervals.
KHASH_MAP_INIT_STR(Region, Region_t*)

// Define our set structure.
typedef struct {
    // Flag to indicate if we wish the set to represent the complement
    //  of genomic intervals.
    bool takeComplement;
    // Our hash table associating chromosomes to intervals.
    khash_t(Region)* regions;
} RegionSet_t;

// A string used to define a set of intervals must follow the following grammar:
// REGIONS -> (REGION,REGIONS) | REGION
// REGION -> STR | STR:INT | STR:INT- | STR:-INT | STR:INT-INT
//  where STR is a string and INT is a non-negative integer.

// Initalize a set.
// Accepts:
//  kstring_t* inputRegions -> Our string that defines intervals. Follows above grammar.
//  bool takeComplement -> Flag to indicate if set should be treated as complement of intervals.
// Returns:
//  RegionSet_t*, The created set.
RegionSet_t* init_region_set(kstring_t* inputRegions, bool takeComplement);

// Used to test if a locus is in the defined set.
// Accepts:
//  RegionSet_t* set -> The set to evaluate membership.
//  kstring_t* chrom -> The chromosome the locus is located on.
//  unsigned int locus -> The coordinate of the locus.
// Returns: bool, If the locus is in the set.
bool query_locus(RegionSet_t* set, kstring_t* chrom, unsigned int locus);

// Used to test if a interval overlaps members of the set.
// Accepts:
//  RegionSet_t* set -> The set to evaluate membership.
//  kstring_t* chrom -> The chromosome the locus is located on.
//  unsigned int startLocus -> The start coordinate of the window.
//  unsigned int endLocus -> The end coordinate of the window.
// Returns: bool, If the window overlaps with set elements.
bool query_overlap(RegionSet_t* set, kstring_t* chrom, unsigned int startLocus, unsigned int endLocus);

// Destroy a set.
// Accepts:
//  RegionSet_t* set -> The set to be destroyed.
// Returns: void.
void destroy_region_set(RegionSet_t* set);

#endif