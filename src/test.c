
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

    int WINDOW_SIZE = 40000, HAP_SIZE = 1, OFFSET_SIZE = 1;

    VCFGenotypeParser* parser = init_vcf_genotype_parser("./data/test.data.40k.vcf");
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