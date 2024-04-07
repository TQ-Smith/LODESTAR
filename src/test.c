
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include <stdio.h>

void print_window_info(Window* window) {
    printf("Window Number: %d\n", window -> winNum);
    printf("Chromosome: %s\n", ks_str(window -> chromosome));
    printf("Window Number on Chromosome: %d\n", window -> winNumOnChrom);
    printf("Start Position: %d\n", window -> startLocus);
    printf("End Position: %d\n", window -> endLocus);
    printf("Number of Loci: %d\n", window -> numLoci);
    printf("\n");
}

int main() {

    int NUM_THREADS = 1;
    int HAP_SIZE = 1, STEP_SIZE = 1, WINDOW_SIZE = 10;

    VCFLocusParser* parser = init_vcf_locus_parser("./data/sliding_window_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> numSamples);

    int numWindows;
    Window** windows = sliding_window(parser, encoder, HAP_SIZE, STEP_SIZE, WINDOW_SIZE, NUM_THREADS, &numWindows);

    for (int i = 0; i < numWindows; i++)
        print_window_info(windows[i]);

    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);
    for (int i = 0; i < numWindows; i++)
        destroy_window(windows[i]);
    free(windows);

    return 0;

}