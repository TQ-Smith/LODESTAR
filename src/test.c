
#include "VCFGenotypeParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include <stdio.h>

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

    int NUM_THREADS = 5;

    int HAP_SIZE = 1, STEP_SIZE = 1, WINDOW_SIZE = 1;

    VCFGenotypeParser* parser = init_vcf_genotype_parser("./data/haplotype_tree_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> num_samples);
    
    int numWindows = 0;
    Window** windows = window_genome(parser, encoder, HAP_SIZE, STEP_SIZE, WINDOW_SIZE, NUM_THREADS, &numWindows);

    printf("\n");
    for (int i = 0; i < numWindows; i++)
        print_window_info(windows[i]);
    
    destroy_vcf_genotype_parser(parser);
    destroy_haplotype_encoder(encoder);
    for (int i = 0; i < numWindows; i++)
        destroy_window(windows[i]);
    free(windows);

    return 0;

}