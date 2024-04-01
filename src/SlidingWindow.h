
// File: SlidingWindow.h
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Slides a window along the contents of a VCF file.

#ifndef _SLIDING_WINDOW_H_
#define _SLIDING_WINDOW_H_

#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "Window.h"

Window** window_genome(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows);

Window* global(VCFLocusParser* parser, HaplotypeEncoder* encoder, int NUM_THREADS);

#endif