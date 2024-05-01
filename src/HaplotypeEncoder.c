
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

    encoder -> labelMap = kh_init(haplotype);

    encoder -> numLeaves = 1;

    return encoder;

}

void relabel_haplotypes(HaplotypeEncoder* encoder) {

    khint_t j;
    khiter_t k;
    int ret;

    for (j = kh_begin(encoder -> labelMap); j != kh_end(encoder -> labelMap); j++) {
		if (!kh_exist(encoder -> labelMap, j)) 
            continue;
		kh_val(encoder -> labelMap, j) = MISSING;
	}
    
    Haplotype newLabel = 0;

    for (int i = 0; i < encoder -> numSamples; i++) {
        if (encoder -> genotypes[i].left == MISSING)
            continue;
        
        k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left);
        if (k == kh_end(encoder -> labelMap)) {
            k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[i].left, &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        if (kh_value(encoder -> labelMap, k) == MISSING) {
            kh_value(encoder -> labelMap, k) = newLabel++;
        }

        encoder -> genotypes[i].left = kh_value(encoder -> labelMap, kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left));

        k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].right);
        if (k == kh_end(encoder -> labelMap)) {
            k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[i].right, &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        if (kh_value(encoder -> labelMap, k) == MISSING) {
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        
        encoder -> genotypes[i].right = kh_value(encoder -> labelMap, kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].right));
    }

    encoder -> numLeaves = newLabel;

}

void add_locus(HaplotypeEncoder* encoder, int numAlleles) {

    for (int i = 0; i < encoder -> numSamples; i++) {
        if (encoder -> numLeaves == 1) {
            encoder -> genotypes[i].left = LEFT_ALLELE(encoder -> locus[i]);
            encoder -> genotypes[i].right = RIGHT_ALLELE(encoder -> locus[i]);
            if (encoder -> genotypes[i].left == numAlleles || encoder -> genotypes[i].right == numAlleles) {
                encoder -> genotypes[i].left = MISSING;
                encoder -> genotypes[i].right = MISSING;
            }
        } else if (encoder -> genotypes[i].left == MISSING || encoder -> genotypes[i].right == MISSING || LEFT_ALLELE(encoder -> locus[i]) == numAlleles || RIGHT_ALLELE(encoder -> locus[i]) == numAlleles) {
            encoder -> genotypes[i].left = MISSING;
            encoder -> genotypes[i].right = MISSING;
        } else {
            encoder -> genotypes[i].left = encoder -> genotypes[i].left * numAlleles + LEFT_ALLELE(encoder -> locus[i]);
            encoder -> genotypes[i].right = encoder -> genotypes[i].right * numAlleles + RIGHT_ALLELE(encoder -> locus[i]);
        }
    }
    
    encoder -> numLeaves = (encoder -> numLeaves) * numAlleles;

    if ((2 * encoder -> numLeaves) / 2 != encoder -> numLeaves)
        relabel_haplotypes(encoder);

}

bool get_next_haplotype(VCFLocusParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE) {

    if (parser -> isEOF)
        return false;
    
    encoder -> chrom -> l = 0;
    ks_overwrite(ks_str(parser -> nextChrom), encoder -> chrom);

    encoder -> startLocus = parser -> nextPos;

    encoder -> numLeaves = 1;

    encoder -> numLoci = 0;

    bool isSameChrom = true;

    int numAlleles;

    while(!(parser -> isEOF) && (encoder -> numLoci < HAP_SIZE) && isSameChrom) {
        get_next_locus(parser, encoder -> chrom, &(encoder -> endLocus), &numAlleles, &(encoder -> locus));
        add_locus(encoder, numAlleles);
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
    kh_destroy(haplotype, encoder -> labelMap);
    free(encoder);
}

// Used to test the haplotype encoder.
/*
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
    get_next_haplotype(parser, encoder, 10);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 10);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);

    printf("\n");
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz");
    printf("\nTest 2\n");
    printf("---------\n");
    printf("Read in 2-loci haplotypes:\n\n");
    get_next_haplotype(parser, encoder, 2);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 2);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 2);
    printf("\nThird Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);

    printf("\n");
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz");
    printf("\nTest 3\n");
    printf("---------\n");
    printf("Read in 1-locus haplotypes:\n\n");
    get_next_haplotype(parser, encoder, 1);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 1);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 1);
    printf("\nThird Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 1);
    printf("\nFourth Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 1);
    printf("\nFifth Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);

    printf("\n");
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test2.vcf.gz");
    printf("\nTest 4\n");
    printf("---------\n");
    get_next_haplotype(parser, encoder, 4);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    printf("\nRelabeled:\n");
    relabel_haplotypes(encoder);
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 10);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    printf("\nRelabeled:\n");
    relabel_haplotypes(encoder);
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);


    printf("\n");

    destroy_haplotype_encoder(encoder);

}
*/