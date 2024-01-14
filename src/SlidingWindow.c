
#include "SlidingWindow.h"

Window* init_window() {
    Window* window = (Window*) calloc(1, sizeof(Window));
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    window -> numLoci = 0;
    return window;
}

Window* get_next_window(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, Window* currentWindow, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    if (parser -> isEOF)
        return NULL;

    Window* nextWindow = init_window();

    kputs(ks_str(parser -> nextChromosome), currentWindow -> chromosome);

    bool isSameChromosome = true;
    int numHaps = currentWindow -> numLoci / HAP_SIZE;

    while(numHaps < WINDOW_SIZE && isSameChromosome) {
        isSameChromosome = get_next_haplotype(parser, encoder, HAP_SIZE);

        // Process haplotype.
        printf("%d\n", encoder -> startLocus);

        currentWindow -> numLoci += encoder -> numLoci;
        numHaps++;
    }
    currentWindow -> endLocus = encoder -> endLocus;
    if (isSameChromosome)
        nextWindow -> numLoci = currentWindow -> numLoci - (HAP_SIZE * OFFSET_SIZE);
    else
        nextWindow -> numLoci = 0;

    return nextWindow;

}

klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    if (parser -> isEOF)
        return NULL;

    klist_t(WindowPtr)* windows = kl_init(WindowPtr);
    
    Window* currentWindow = init_window();
    Window* nextWindow = NULL;
    Window* temp = NULL;

    while ((nextWindow = get_next_window(parser, encoder, currentWindow, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE)) != NULL) {

        // Process currentWindow.

        *kl_pushp(WindowPtr, windows) = currentWindow;

        temp = currentWindow;
        currentWindow = nextWindow;
        nextWindow = temp;

    }

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
    printf("Chromosome: %s\n", ks_str(window -> chromosome));
    printf("Start Position: %d\n", window -> startLocus);
    printf("End Position: %d\n", window -> endLocus);
    printf("Number of Loci: %d\n", window -> numLoci);
    printf("\n");
}

int main() {

    int WINDOW_SIZE = 2, HAP_SIZE = 1, OFFSET_SIZE = 1;

    VCFGenotypeParser* parser = init_vcf_genotype_parser("sliding_window_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> num_samples);
    
    klist_t(WindowPtr)* windows = slide_through_genome(parser, encoder, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE);
    
    printf("\nHaplotype Size of %d SNPs\nOffset Size of %d Haplotypes\nWindow Size of %d Haplotypes\n", HAP_SIZE, OFFSET_SIZE, WINDOW_SIZE);

    printf("\nWindows:");
    printf("\n-------\n\n");
    int i = 1;
    for (kliter_t(WindowPtr)* it = kl_begin(windows); it != kl_end(windows); it = kl_next(it)) {
        printf("Window %d:\n", i);
        print_window_info(kl_val(it));
        i++;
    }

    kl_destroy(WindowPtr, windows);
    destroy_vcf_genotype_parser(parser);
    destroy_haplotype_encoder(encoder);

    return 0;

}