
#include "SlidingWindow.h"

Window* init_window() {
    Window* window = (Window*) calloc(1, sizeof(Window));
    window -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    window -> numLoci = 0;
    return window;
}

Window* get_next_window(VCFGenotypeParser* parser, HaplotypeTree* tree, Window* currentWindow, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    if (parser -> isEOF)
        return NULL;

    Window* nextWindow = init_window();

    bool noNextHapOnSameChromosome = true;
    int numHaps = currentWindow -> numLoci / HAP_SIZE;

    while(noNextHapOnSameChromosome && numHaps < WINDOW_SIZE) {
        noNextHapOnSameChromosome = get_next_haplotype(parser, tree, HAP_SIZE);

        printf("Hap\n");

        currentWindow -> numLoci += tree -> numLoci;
        if (numHaps % OFFSET_SIZE == 0)
            nextWindow -> startLocus = parser -> nextPosition;
        numHaps++;
    }
    currentWindow -> endLocus = tree -> endLocus;
    if (strcmp(ks_str(currentWindow -> chromosome), ks_str(parser -> nextChromosome)) == 0)
        nextWindow -> numLoci = 0;
    else
        nextWindow -> numLoci = currentWindow -> numLoci - (HAP_SIZE * OFFSET_SIZE);
    kputs(ks_str(parser -> nextChromosome), nextWindow -> chromosome);

    return nextWindow;

}

klist_t(WindowPtr)* slide_through_genome(VCFGenotypeParser* parser, HaplotypeTree* tree, int WINDOW_SIZE, int HAP_SIZE, int OFFSET_SIZE) {
    
    if (parser -> isEOF)
        return NULL;

    klist_t(WindowPtr)* windows = kl_init(WindowPtr);
    
    Window* currentWindow = init_window();
    Window* nextWindow = NULL;
    Window* temp = NULL;

    currentWindow -> startLocus = parser -> nextPosition;
    kputs(ks_str(parser -> nextChromosome), currentWindow -> chromosome);
    currentWindow -> numLoci = 0;
    
    while ((nextWindow = get_next_window(parser, tree, currentWindow, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE)) != NULL) {

        // Process current Window.

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

    int WINDOW_SIZE = 2, HAP_SIZE = 2, OFFSET_SIZE = 1;

    VCFGenotypeParser* parser = init_vcf_genotype_parser("haplotype_tree_test.vcf.gz");
    HaplotypeTree* tree = init_haplotype_tree(parser -> num_samples);
    
    klist_t(WindowPtr)* windows = slide_through_genome(parser, tree, WINDOW_SIZE, HAP_SIZE, OFFSET_SIZE);
    
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
    destroy_haplotype_tree(tree);

    return 0;

}