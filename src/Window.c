
// File: Window.c
// Date: 8 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines the attributes of a window in the genome.

#include "Window.h"
#include <stdlib.h>

Window_t* init_window() {
    // Allocate memory and set defaults.
    Window_t* window = (Window_t*) calloc(1, sizeof(Window_t));
    window -> chromosome = NULL;
    window -> X = NULL;
    window -> x0 = NULL;
    window -> saveIBS = false;
    window -> ibs = NULL;
    // window -> pval = -1;
    window -> t = -1;
    window -> sigma = -1;
    window -> dropWindow = false;
    return window;
}

void destroy_window(Window_t* window, int n) {
    // Free all memory when not NULL.
    if (window == NULL)
        return;
    if (window -> X != NULL) {
        for (int i = 0; i < n; i++)
            free(window -> X[i]);
        free(window -> X);
    }
    if (window -> x0 != NULL)
        free(window -> x0);
    if (window -> ibs != NULL)
        free(window -> ibs);
    destroy_kstring(window -> chromosome);
    free(window);
    window = NULL;
}