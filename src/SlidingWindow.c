
#include "SlidingWindow.h"

Window* init_window() {
    Window* window = (Window*) calloc(1, sizeof(Window));
    window -> windowNum = 1;
    window -> windowNumOnChromosome = 1;
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    window -> numLoci = 0;
    return window;
}

Window* get_next_window(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, Window* currentWindow, int* startLoci, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    if (parser -> isEOF)
        return NULL;

    Window* nextWindow = init_window();

    kputs(ks_str(parser -> nextChromosome), currentWindow -> chromosome);

    bool isSameChromosome = true;
    int numHapsInOverlap = currentWindow -> numLoci / HAP_SIZE;

    while(numHapsInOverlap < WINDOW_SIZE && isSameChromosome) {
        isSameChromosome = get_next_haplotype(parser, encoder, HAP_SIZE);
        // Process haplotype.
        if (numHapsInOverlap % OFFSET_SIZE == 0)
            startLoci[(currentWindow -> windowNumOnChromosome + numHapsInOverlap) % ((WINDOW_SIZE - OFFSET_SIZE) / OFFSET_SIZE + 1)] = encoder -> startLocus;
        currentWindow -> numLoci += encoder -> numLoci;
        numHapsInOverlap++;
    }
    currentWindow -> startLocus = startLoci[currentWindow -> windowNumOnChromosome % ((WINDOW_SIZE - OFFSET_SIZE) / OFFSET_SIZE + 1)];
    currentWindow -> endLocus = encoder -> endLocus;
    nextWindow -> windowNum = currentWindow -> windowNum + 1;
    if (isSameChromosome) {
        nextWindow -> numLoci = currentWindow -> numLoci - (HAP_SIZE * OFFSET_SIZE);
        nextWindow -> windowNumOnChromosome = currentWindow -> windowNumOnChromosome + 1;
    } else {
        nextWindow -> numLoci = 0;
        nextWindow -> windowNumOnChromosome = 1;
    }

    return nextWindow;

}

klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    if (parser -> isEOF)
        return NULL;

    klist_t(WindowPtr)* windows = kl_init(WindowPtr);

    int* startLoci = (int*) calloc((WINDOW_SIZE - OFFSET_SIZE) / OFFSET_SIZE + 1, sizeof(int));
    
    Window* currentWindow = init_window();
    Window* nextWindow = NULL;
    Window* temp = NULL;

    int i = 1;
    while ((nextWindow = get_next_window(parser, encoder, currentWindow, startLoci, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE)) != NULL) {

        // Process currentWindow.
        printf("Processed Window %d!\n", i++);

        *kl_pushp(WindowPtr, windows) = currentWindow;

        temp = currentWindow;
        currentWindow = nextWindow;
        nextWindow = temp;

    }

    destroy_window(currentWindow);
    free(startLoci);

    return windows;

}

void destroy_window(Window* window) {
    if (window == NULL)
        return;
    free(ks_str(window -> chromosome)); 
    free(window -> chromosome);
    free(window);
}

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