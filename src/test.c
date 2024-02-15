
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

    VCFGenotypeParser* parser = init_vcf_genotype_parser("./data/haplotype_tree_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> num_samples);
    
    
    destroy_vcf_genotype_parser(parser);
    destroy_haplotype_encoder(encoder);

    return 0;

}