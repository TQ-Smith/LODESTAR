
#include "Window.h"

#include <stdlib.h>

Window* init_window() {
    // Allocate the structure.
    Window* window = (Window*) calloc(1, sizeof(Window));
    // Set default numbering.
    window -> windowNum = 1;
    window -> windowNumOnChromosome = 1;
    window -> numLoci = 0;
    // Allocate string for chromosome.
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    // Return window.
    return window;
}

void destroy_window(Window* window) {
    // Cannot destroy a NULL window.
    if (window == NULL)
        return;

    // Free the chromosome string.
    free(ks_str(window -> chromosome)); free(window -> chromosome);

    // Free the structure.
    free(window);
}