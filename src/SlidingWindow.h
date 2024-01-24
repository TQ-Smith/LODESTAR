
// File: SlidingWindow.h
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Slides a window along the contents of a VCF file.

#ifndef _SLIDING_WINDOW_
#define _SLIDING_WINDOW_

#include "VCFGenotypeParser.h"

#include "HaplotypeEncoder.h"

#include "Window.h"

#include "ASD.h"

// We are using klib list.
#include "../klib/klist.h"

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
klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, ASD* asd, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE);

#endif