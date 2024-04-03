
// File: HaplotypeEncoder.c
// Date: 
// Author: TQ Smith
// Purpose: 

#include "HaplotypeEncoder.h"

#include <math.h>

HaplotypeEncoder* init_haplotype_encoder(int numSamples) {

    HaplotypeEncoder* encoder = (HaplotypeEncoder*) calloc(1, sizeof(HaplotypeEncoder));

    encoder -> numSamples = numSamples;
    encoder -> locus = (Locus*) calloc(numSamples, sizeof(Locus));
    encoder -> genotypes = (Genotype*) calloc(numSamples, sizeof(Genotype));

    encoder -> chrom = (kstring_t*) calloc(1, sizeof(kstring_t));

    encoder -> numLeaves = 1;

    encoder -> labelMap = kh_init(int64);

    return encoder;

}

void relabel_haplotypes(HaplotypeEncoder* encoder) {

    khint_t k;
    int ret, newLabel = 0;
    
    for (k = kh_begin(encoder -> labelMap); k != kh_end(encoder -> labelMap); k++)
        kh_value(encoder -> labelMap, k) = MISSING;

    for (int i = 0; i < encoder -> numSamples; i++) {

    }

    encoder -> numLeaves = newLabel;

}

void add_locus(HaplotypeEncoder* encoder, int numAlleles, bool collapseMissingGenotypes) {

    //if (encoder -> numLeaves * numAlleles == MISSING || (encoder -> numLeaves * numAlleles) / numAlleles != encoder -> numLeaves)
        //relabel_haplotypes(encoder);
    
    if (encoder -> numLeaves >= 2)
        relabel_haplotypes(encoder);

    for (int i = 0; i < encoder -> numSamples; i++) {
        if (encoder -> numLeaves == 1) {
            encoder -> genotypes[i].left = LEFT_ALLELE(encoder -> locus[i]);
            encoder -> genotypes[i].right = RIGHT_ALLELE(encoder -> locus[i]);
            if (collapseMissingGenotypes && (encoder -> genotypes[i].left == numAlleles || encoder -> genotypes[i].right == numAlleles)) {
                encoder -> genotypes[i].left = MISSING;
                encoder -> genotypes[i].right = MISSING;
            }
        } else if (collapseMissingGenotypes && (encoder -> genotypes[i].left == MISSING || encoder -> genotypes[i].right == MISSING || LEFT_ALLELE(encoder -> locus[i]) == numAlleles || RIGHT_ALLELE(encoder -> locus[i]) == numAlleles)) {
            encoder -> genotypes[i].left = MISSING;
            encoder -> genotypes[i].right = MISSING;
        } else {
            encoder -> genotypes[i].left = encoder -> genotypes[i].left * numAlleles + LEFT_ALLELE(encoder -> locus[i]);
            encoder -> genotypes[i].right = encoder -> genotypes[i].right * numAlleles + RIGHT_ALLELE(encoder -> locus[i]);
        }
    }
    
    encoder -> numLeaves = (encoder -> numLeaves) * numAlleles;

}

bool get_next_haplotype(VCFLocusParser* parser, HaplotypeEncoder* encoder, bool collapseMissingGenotypes, int HAP_SIZE) {

    if (parser -> isEOF)
        return false;
    
    encoder -> chrom -> l = 0;
    kputs(ks_str(parser -> nextChrom), encoder -> chrom);

    encoder -> startLocus = parser -> nextPos;

    encoder -> numLeaves = 1;

    encoder -> numLoci = 0;

    bool isSameChrom = true;

    int numAlleles;

    while(!(parser -> isEOF) && (encoder -> numLoci < HAP_SIZE) && isSameChrom) {
        get_next_locus(parser, encoder -> chrom, &(encoder -> endLocus), &numAlleles, &(encoder -> locus));
        add_locus(encoder, numAlleles, collapseMissingGenotypes);
        isSameChrom = strcmp(ks_str(encoder -> chrom), ks_str(parser -> nextChrom)) == 0;
        encoder -> numLoci++;
    }

    return !(parser -> isEOF) && encoder -> numLoci == HAP_SIZE && isSameChrom;

}

void destroy_haplotype_encoder(HaplotypeEncoder* encoder) {
    if (encoder == NULL)
        return;
    free(encoder -> locus);
    free(encoder -> genotypes);
    free(ks_str(encoder -> chrom)); free(encoder -> chrom);
    kh_destroy(int64, encoder -> labelMap);
    free(encoder);
}

// Used to test the haplotype encoder.


void print_encoder_info(HaplotypeEncoder* encoder) {
    printf("Chromosome: %s\n", ks_str(encoder -> chrom));
    printf("Start locus: %d\n", encoder -> startLocus);
    printf("End locus: %d\n", encoder -> endLocus);
    printf("Number of loci: %d\n", encoder -> numLoci);
    printf("Sample Haplotypes:\n");
    for (int i = 0; i < encoder -> numSamples; i++)
        printf("Sample %d -> %ld/%ld\n", i + 1, encoder -> genotypes[i].left, encoder -> genotypes[i].right);
}

int main() {

    VCFLocusParser* parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> numSamples);
    printf("\nTest 1\n");
    printf("---------\n");
    printf("Read in whole chromosomes:\n\n");
    get_next_haplotype(parser, encoder, true, 10);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 10);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);

    printf("\n");
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz");
    printf("\nTest 2\n");
    printf("---------\n");
    printf("Read in 2-loci haplotypes:\n\n");
    get_next_haplotype(parser, encoder, true, 2);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 2);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 2);
    printf("\nThird Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);

    printf("\n");
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz");
    printf("\nTest 3\n");
    printf("---------\n");
    printf("Read in 1-locus haplotypes:\n\n");
    get_next_haplotype(parser, encoder, true, 1);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 1);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 1);
    printf("\nThird Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 1);
    printf("\nFourth Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, true, 1);
    printf("\nFifth Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);

    printf("\n");

    destroy_haplotype_encoder(encoder);

}
