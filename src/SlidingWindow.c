
// File: SlidingWindow.c
// Date: 
// Author: TQ Smith
// Purpose: 

#include "SlidingWindow.h"

#include "Matrix.h"

#include "../klib/klist.h"
#define do_not_free_window_ptr(w)
KLIST_INIT(WindowPtr, Window*, do_not_free_window_ptr)

#include <pthread.h>

#define SWAP(a, b, TEMP) TEMP = a; a = b; b = TEMP

pthread_mutex_t vcfLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t windowListLock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t globalLock = PTHREAD_MUTEX_INITIALIZER;

Window* get_next_window() {

}

void* sliding_window_multi_thread(void* arg) {

}

void sliding_window_single_thread() {

}

Window** window_genome(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE, int STEP_SIZE, int WINDOW_SIZE, int NUM_THREADS, int* numWindows) {

}

Window* global(VCFLocusParser* parser, HaplotypeEncoder* encoder, int NUM_THREADS) {

}