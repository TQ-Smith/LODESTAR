
// File: SlidingWindow.h
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Slides a window along the contents of a VCF file.

#ifndef _SLIDING_WINDOW_
#define _SLIDING_WINDOW_

#include "VCFGenotypeParser.h"

#include "HaplotypeEncoder.h"

#include "Window.h"

Window** window_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows);

#endif