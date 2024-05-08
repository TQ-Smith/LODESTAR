
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
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    window -> X = NULL;
    window -> saveIBS = false;
    window -> ibs = NULL;
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
    if (window -> ibs != NULL)
        free(window -> ibs);
    free(ks_str(window -> chromosome)); free(window -> chromosome);
    free(window);
}