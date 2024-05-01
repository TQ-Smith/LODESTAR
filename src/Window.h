
#ifndef _WINDOW_H_
#define _WINDOW_H_

#include "../lib/kstring.h"

#include "AlleleSharingDistance.h"

typedef struct {

    int winNum;
    int winNumOnChrom;
    kstring_t* chromosome;
    unsigned int startLocus;
    unsigned int endLocus;
    int numLoci;

    double** X;
    IBS* ibs;

} Window;

Window* init_window();

void destroy_window(Window* window, int n);

#endif