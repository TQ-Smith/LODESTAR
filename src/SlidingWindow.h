
// File: SlidingWindow.h
// Date: 8 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Perform ASD and MDS along windows in the genome or just genome-wide.

#ifndef _SLIDING_WINDOW_H_
#define _SLIDING_WINDOW_H_

#include "HaplotypeEncoder.h"
#include "Window.h"

// Perform ASD and MDS along windows in the genome.
// Accepts:
//  VCFLocusParser_t* parser -> The configured parser to read in a VCF.
//  HaplotypeEncoder_t* encoder -> Used to encode haplotypes while reading in the VCF.
//  int k -> The dimension to project samples into after calculating ASD matrix.
//  int HAP_SIZE -> The maximum number of loci that can compose a haplotype.
//  int STEP_SIZE -> The number of haplotypes to advance the start of the current
//                      window to determine the start of the next window. Must be
//                      less than or equal to WINDOW_SIZE.
//  int WINDOW_SIZE -> The maximum number of haplotypes that can compose a window.
//  int NUM_THREADS -> The number of threads to use in the calculations.
//  int* numWindows -> Sets the total number of processed windows, including the genome-wide.
// Returns: Window_t**, An array of window pointers.
Window_t** sliding_window(VCFLocusParser_t* parser, HaplotypeEncoder_t* encoder, int k, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows);

// Perform ASD and MDS just for the genome-wide window.
// Accepts:
//  VCFLocusParser_t* parser -> The configured parser to read in a VCF.
//  HaplotypeEncoder_t* encoder -> Used to encode haplotypes while reading in the VCF.
//  int k -> The dimension to project samples into after calculating ASD matrix.
//  int HAP_SIZE -> The maximum number of loci that can compose a haplotype.
//  int NUM_THREADS -> The number of threads to use in the calculations.
// Returns: Window_t*, The genome-wide window.
Window_t* global_window(VCFLocusParser_t* parser, HaplotypeEncoder_t* encoder, int k, int HAP_SIZE, int NUM_THREADS);

#endif