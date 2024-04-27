
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include "ProcrustesAnalysis.h"

#include <stdio.h>

#include <stdlib.h>

void print_window_info(Window* window, int n, int k) {
    printf("Window Number: %d\n", window -> winNum);
    printf("Chromosome: %s\n", ks_str(window -> chromosome));
    printf("Window Number on Chromosome: %d\n", window -> winNumOnChrom);
    printf("Start Position: %d\n", window -> startLocus);
    printf("End Position: %d\n", window -> endLocus);
    printf("Number of Loci: %d\n", window -> numLoci);
    printf("X = \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++)
            printf("%5lf\t", window -> X[i][j]);
        printf("\n");
    }
    printf("\n\n");
}

int main () {

    int k = 1;
    int NUM_THREADS = 1;
    int HAP_SIZE = 1, STEP_SIZE = 3, WINDOW_SIZE = 4;

    VCFLocusParser* parser = init_vcf_locus_parser("./data/sliding_window_test2.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> numSamples);

    int numWindows;
    Window** windows = sliding_window(parser, encoder, k, HAP_SIZE, STEP_SIZE, WINDOW_SIZE, NUM_THREADS, &numWindows);

    for (int i = 0; i < numWindows; i++)
        print_window_info(windows[i], encoder -> numSamples, k);

    for (int i = 0; i < numWindows; i++)
        destroy_window(windows[i], parser -> numSamples);
    free(windows);
    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);

    return 0;

}