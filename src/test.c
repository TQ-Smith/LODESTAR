
#include "VCFGenotypeParser.h"

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

    int HAP_SIZE = 1, STEP_SIZE = 1, WINDOW_SIZE = 3;

    VCFGenotypeParser* parser = init_vcf_genotype_parser("./data/sliding_window_test2.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> numSamples);

    return 0;

}