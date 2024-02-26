
#ifndef _WINDOW_
#define _WINDOW_

#include "../klib/kstring.h"

// A structure to hold a window's information.
//  Will change significantly given the application.
typedef struct {

    // The window number out of all processed windows.
    int winNum;
    // The window number on a griven chromosome.
    int winNumOnChrom;
    // The chromosome the window is on.
    kstring_t* chromosome;
    // The start locus of the window.
    int startLocus;
    // The end locus of the window.
    int endLocus;
    // The number of loci within the window.
    int numLoci;

} Window;

// Creates a window object.
//  Will change with the given application.
// Accepts:
//  void.
// Returns:
//  Window*, A pointer to a new window structure.
Window* init_window();

// Deallocates the memory occupied by a window.
//  Will change with the given application.
// Accepts:
//  Window* window -> The window to deallocate.
// Returns:
//  void.
void destroy_window(Window* window);

#endif