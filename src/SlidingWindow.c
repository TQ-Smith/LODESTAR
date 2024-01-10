
#include "SlidingWindow.h"

Window* init_window() {
    Window* window = (Window*) calloc(1, sizeof(Window));
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    window -> numLoci = 0;
    return window;
}

Window* get_next_window(VCFGenotypeParser* parser, HaplotypeTree* tree, Window* current_window, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    return NULL;
}

klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeTree* tree, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    if (parser -> isEOF)
        return NULL;

    klist_t(WindowPtr)* windows = kl_init(WindowPtr);

    Window* current_window = init_window();
    Window* next_window = NULL;
    Window* temp = NULL;

    while(!(parser -> isEOF)) {
        next_window = get_next_window(parser, tree, current_window, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE);
        // Process current and then save current to list.
        *kl_pushp(WindowPtr, windows) = current_window;

        temp = next_window;
        next_window = current_window;
        current_window = temp;
    }

    return windows;

}

void destroy_window(Window* window) {
    free(ks_str(window -> chromosome)); free(window -> chromosome);
    free(window);
}

int main() {

}