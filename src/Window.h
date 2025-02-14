
// File: Window.h
// Date: 8 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines the attributes of a window in the genome.

#ifndef _WINDOW_H_
#define _WINDOW_H_

#include "../lib/kstring.h"
#include "IdentityByState.h"
#include <stdbool.h>

typedef struct {

    // The window number along the whole chromosome.
    int winNum;
    // The window number on the current chromosome.
    int winNumOnChrom;
    // The chromosome the window is on.
    kstring_t* chromosome;
    // The start coordinate of the first haplotype in the window.
    unsigned int startCoord;
    // The end coordinate of the last haplotype in the window.
    unsigned int endCoord;
    // The number of loci in the window.
    int numLoci;
    // The number of haplotypes in the window.
    int numHaps;
    
    // The matrix of points representing samples in k-dimensional space.
    double** X;
    // The vector of X's column averages.
    double* x0;

    // A flag to indicate if the IBS counts for the window should be saved.
    bool saveIBS;
    // The IBS counts for the window. If saveIBS is not set, pointer is NULL.
    IBS_t* ibs;

    // The p-value associated with the window.
    double pval;

    // The t-statistic associated with the window.
    double t;

    // A flag to indicate if we drop the window because it exceeds gap threshold.
    bool dropWindow;

} Window_t;


// Create and allocate memory for a window.
// Accepts: void.
// Returns: Window_t*, The pointer to the created window.
Window_t* init_window();

// Destroy memory occupied by a window.
// Accepts:
//  Window_t* window -> Pointer to window structure to destroy.
//  int n -> The number of samples in the set of points.
//              Used to free window -> X.
// Returns: void.
void destroy_window(Window_t* window, int n);

#endif