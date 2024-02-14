
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

void print_encoder_info(HaplotypeEncoder* encoder) {
    printf("Chromosome: %s\n", ks_str(encoder -> chromosome));
    printf("Start locus: %d\n", encoder -> startLocus);
    printf("End locus: %d\n", encoder -> endLocus);
    printf("Number of loci: %d\n", encoder -> numLoci);
    printf("Sample Haplotypes:\n");
    for (int i = 0; i < encoder -> numSamples; i++)
        printf("Sample %d -> %lf/%lf\n", i + 1, encoder -> leftHaplotype[i], encoder -> rightHaplotype[i]);
}

int main() {

    VCFGenotypeParser* parser = init_vcf_genotype_parser("./data/haplotype_tree_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> num_samples);
    
    printf("\nTest 1\n");
    printf("---------\n");
    printf("Read in whole chromosomes:\n\n");
    get_next_haplotype(parser, encoder, false, 10);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, false, 10);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    parser = init_vcf_genotype_parser("./data/haplotype_tree_test.vcf.gz");
    printf("\nTest 2\n");
    printf("---------\n");
    printf("Read in 2-loci haplotypes:\n\n");
    get_next_haplotype(parser, encoder, false, 2);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, false, 2);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, false, 2);
    printf("\nThird Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    parser = init_vcf_genotype_parser("./data/haplotype_tree_test.vcf.gz");
    printf("\nTest 3\n");
    printf("---------\n");
    printf("Read in 1-locus haplotypes:\n\n");
    get_next_haplotype(parser, encoder, false, 1);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, false, 1);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, false, 1);
    printf("\nThird Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, false, 1);
    printf("\nFourth Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, false, 1);
    printf("\nFifth Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    parser = init_vcf_genotype_parser("./data/haplotype_tree_test.vcf.gz");
    printf("\nTest 4\n");
    printf("---------\n");
    printf("Read in whole chromosomes and collapse genotypes:\n\n");
    get_next_haplotype(parser, encoder, true, 10);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 10);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    
    destroy_haplotype_encoder(encoder);

    return 0;

}