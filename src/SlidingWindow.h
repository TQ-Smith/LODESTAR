
#ifndef _SLIDING_WINDOW_
#define _SLIDING_WINDOW_

#include "VCFGenotypeParser.h"

#include "HaplotypeEncoder.h"

#include "../klib/klist.h"

typedef struct {

    kstring_t* chromosome;
    int startLocus;
    int endLocus;
    int numLoci;

} Window;

Window* init_window();

void destroy_window(Window* window);

#define destroy_w(w) destroy_window((w) -> data)
KLIST_INIT(WindowPtr, Window*, destroy_w)

klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE);

#endif