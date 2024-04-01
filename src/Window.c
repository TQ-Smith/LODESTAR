
#include "Window.h"

#include <stdlib.h>

Window* init_window() {
    Window* window = (Window*) calloc(1, sizeof(Window));
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    return window;
}

void destroy_window(Window* window) {
    if (window == NULL)
        return;
    free(ks_str(window -> chromosome)); free(window -> chromosome);
    free(window);
}