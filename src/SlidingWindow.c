
// File: SlidingWindow.c
// Date: 18 Janurary 2023
// Author: TQ Smith
// Purpose: Slides a window along the contents of a VCF file.

#include "SlidingWindow.h"

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

// A method to get the next window in sliding process.
// Accepts:
//  VCFGenotypeParser* parser -> The VCF parser.
//  HaplotypeEncoder* encoder -> The encoder used to label haplotypes.
//  Window* currentWindow -> The window currently being processed. Needed for overlap
//                              calculations to crete the next window.
//  int* startLoci -> An array to hold the start loci of the next windows within the current window.
//  int WINDOW_SIZE -> The number of haplotypes in the window.
//  int HAP_SIZE -> The number of loci in a haplotype.
//  int OFFSET_SIZE -> The number of haplotypes in the offset.
Window* get_next_window(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, Window* currentWindow, int* startLoci, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    // If EOF, there is no new window.
    if (parser -> isEOF)
        return NULL;
    
    // Allocate the next window.
    Window* nextWindow = init_window();

    // Current window lies on the chromosome of the next VCF record.
    kputs(ks_str(parser -> nextChromosome), currentWindow -> chromosome);

    // Flag used to determine if next haplotype is on the same chromsome.
    bool isSameChromosome = true;

    // Number of haplotypes within the overlap from the previous window.
    int numHapsInOverlap = currentWindow -> numLoci / HAP_SIZE;

    // Read in the next window.
    while(numHapsInOverlap < WINDOW_SIZE && isSameChromosome) {
        // Get the next haplotype.
        isSameChromosome = get_next_haplotype(parser, encoder, true, HAP_SIZE);

        // Process haplotype.

        // If the haplotype encountered is a start position for a future window, save the haplotype's start position.
        if (numHapsInOverlap % OFFSET_SIZE == 0)
            startLoci[(currentWindow -> windowNumOnChromosome + numHapsInOverlap) % ((WINDOW_SIZE - OFFSET_SIZE) / OFFSET_SIZE + 1)] = encoder -> startLocus;
        currentWindow -> numLoci += encoder -> numLoci;
        numHapsInOverlap++;
    }
    // Set start locus of the currentWindow.
    currentWindow -> startLocus = startLoci[currentWindow -> windowNumOnChromosome % ((WINDOW_SIZE - OFFSET_SIZE) / OFFSET_SIZE + 1)];
    // Set end locus of the currentWindow.
    currentWindow -> endLocus = encoder -> endLocus;
    // Set window number for the next window.
    nextWindow -> windowNum = currentWindow -> windowNum + 1;
    // If the next window is on the same chromsome ...
    if (isSameChromosome) {
        // Subtract number of loci within the offset.
        nextWindow -> numLoci = currentWindow -> numLoci - (HAP_SIZE * OFFSET_SIZE);
        // Next window on the same chromosome.
        nextWindow -> windowNumOnChromosome = currentWindow -> windowNumOnChromosome + 1;
    } else {
        // If next window is on a different chromosome, then reset counters.
        nextWindow -> numLoci = 0;
        nextWindow -> windowNumOnChromosome = 1;
    }

    // Return the newly allocated window.
    return nextWindow;

}

klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    // If EOF, there are no windows to process.
    if (parser -> isEOF)
        return NULL;
    
    // Create our list of window pointers.
    klist_t(WindowPtr)* windows = kl_init(WindowPtr);

    // Allocate our list of start locations.
    int* startLoci = (int*) calloc((WINDOW_SIZE - OFFSET_SIZE) / OFFSET_SIZE + 1, sizeof(int));
    
    // Create the first window.
    Window* currentWindow = init_window();

    // Will hold pointer to next window.
    Window* nextWindow = NULL;

    // Used to swap curretWindow and nextWindow.
    Window* temp = NULL;

    // While there is a window to process.
    while ((nextWindow = get_next_window(parser, encoder, currentWindow, startLoci, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE)) != NULL) {

        // Process currentWindow.
        
        // Add the currentWindow to the list.
        *kl_pushp(WindowPtr, windows) = currentWindow;

        // Swap currentWindow and nextWindow.
        temp = currentWindow;
        currentWindow = nextWindow;
        nextWindow = temp;

    }

    // The last window becomes the currentWindow, is unused, and not added to list.
    //  Free unused window.
    destroy_window(currentWindow);

    // Free startLoci array.
    free(startLoci);

    // Return the structure.
    return windows;

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

// Used to test sliding window.

void print_window_info(Window* window) {
    printf("Window Number: %d\n", window -> windowNum);
    printf("Chromosome: %s\n", ks_str(window -> chromosome));
    printf("Window Number on Chromosome: %d\n", window -> windowNumOnChromosome);
    printf("Start Position: %d\n", window -> startLocus);
    printf("End Position: %d\n", window -> endLocus);
    printf("Number of Loci: %d\n", window -> numLoci);
    printf("\n");
}

int main() {

    int WINDOW_SIZE = 10, HAP_SIZE = 100, OFFSET_SIZE = 1;

    VCFGenotypeParser* parser = init_vcf_genotype_parser("sliding_window_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> num_samples);
    
    klist_t(WindowPtr)* windows = slide_through_genome(parser, encoder, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE);
    
    printf("\nHaplotype Size of %d SNPs\nOffset Size of %d Haplotypes\nWindow Size of %d Haplotypes\n", HAP_SIZE, OFFSET_SIZE, WINDOW_SIZE);

    printf("\nWindows:");
    printf("\n-------\n\n");
    for (kliter_t(WindowPtr)* it = kl_begin(windows); it != kl_end(windows); it = kl_next(it)) {
        print_window_info(kl_val(it));
    }

    kl_destroy(WindowPtr, windows);
    destroy_vcf_genotype_parser(parser);
    destroy_haplotype_encoder(encoder);

    return 0;

}