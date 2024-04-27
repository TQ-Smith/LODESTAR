
#include "Window.h"

#include <stdlib.h>

Window* init_window() {
    Window* window = (Window*) calloc(1, sizeof(Window));
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    window -> X = NULL;
    window -> ibs = NULL;
    return window;
}

void destroy_window(Window* window, int n) {
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