
// File: SlidingWindow.h
// Date: 18 Janurary 2023
// Author: TQ Smith
// Purpose: Slides a window along the contents of a VCF file.

#ifndef _SLIDING_WINDOW_
#define _SLIDING_WINDOW_

#include "VCFGenotypeParser.h"

#include "HaplotypeEncoder.h"

// We are using klib list.
#include "../klib/klist.h"

// A structure to hold a window's information.
//  Will change significantly given the application.
typedef struct {

    // The window number out of all processed windows.
    int windowNum;
    // The window number on a griven chromosome.
    int windowNumOnChromosome;
    // The chromosome the window is on.
    kstring_t* chromosome;
    // The start locus of the window.
    int startLocus;
    // The end locus of the window.
    int endLocus;
    // The number of loci within the window.
    int numLoci;

} Window;

// Creates a window object.
//  Will change with the given application.
// Accepts:
//  void.
// Returns:
//  Window*, A pointer to a new window structure.
Window* init_window();

// Deallocates the memory occupied by a window.
//  Will change with the given application.
// Accepts:
//  Window* window -> The window to deallocate.
// Returns:
//  void.
void destroy_window(Window* window);

// A macro is used to wrap the destroy_window function when
//  deallocating a list of windows. Will change with the 
//  given application.
#define destroy_w(w) destroy_window((w) -> data)
KLIST_INIT(WindowPtr, Window*, destroy_w)

// Method to slide through window and generate a list of windows.
// Accepts:
//  VCFGenotpyeParser* parser -> The VCF file parser to read.
//  HaplotypeEncoder* encoder -> The encoder used to encode haplotypes.
//  int WINDOW_SIZE -> The number of haplotypes in a window.
//  int HAP_SIZE -> The number of loci in a haplotype.
//  int OFFSET_SIZE -> The number of haplotypes in the offset.
// Returns:
//  klist_t(WindowPtr)*, A pointer to a klist of window pointers.
klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE);

#endif