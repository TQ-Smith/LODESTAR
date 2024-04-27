
// File: 
// Date: 
// Author: 
// Purpose: 

#ifndef _SLIDING_WINDOW_H_
#define _SLIDING_WINDOW_H_

#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "Window.h"

Window** sliding_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int k, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows);

Window* global_window(VCFLocusParser* parser, HaplotypeEncoder* encoder, int k, int HAP_SIZE, int NUM_THREADS);

#endif