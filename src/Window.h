
#ifndef _WINDOW_H_
#define _WINDOW_H_

#include "../klib/kstring.h"

typedef struct {

    int winNum;
    int winNumOnChrom;
    kstring_t* chromosome;
    int startLocus;
    int endLocus;
    int numLoci;

} Window;

Window* init_window();

void destroy_window(Window* window);

#endif