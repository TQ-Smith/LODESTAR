
// File: HaplotypeEncoder.c
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Track the haplotypes of each sample using Godel Encoding.

#include "HaplotypeEncoder.h"

#include <math.h>

// Macros to get the left and right allele from the 1-byte encoding.
#define LEFT_ALLELE(a) (a >> 4)
#define RIGHT_ALLELE(a) (a & 0x0F)

#define MISSING (NAN)

HaplotypeEncoder* init_haplotype_encoder(int numSamples) {

    // Allocate the structure's memory.
    HaplotypeEncoder* encoder = (HaplotypeEncoder*) calloc(1, sizeof(HaplotypeEncoder));

    // Allocate memory used for arrays.
    encoder -> numSamples = numSamples;
    encoder -> genos = (GENOTYPE*) calloc(numSamples, sizeof(GENOTYPE));
    encoder -> leftHaps = (double*) calloc(numSamples, sizeof(double));
    encoder -> rightHaps = (double*) calloc(numSamples, sizeof(double));

    // Allocate string to hold chromosome.
    encoder -> chrom = (kstring_t*) calloc(1, sizeof(kstring_t));

    for (int i = 0; i < MAX_NUMBER_OF_LOCI; i++) {
        LOG_PRIMES[i] = log(LOG_PRIMES[i]);
    }

    // Return the encoder.
    return encoder;

}

void add_locus(HaplotypeEncoder* encoder, int numAlleles, bool collapseMissingGenotypes) {

    for (int i = 0; i < encoder -> numSamples; i++) {
        if (encoder -> numLoci == 0) {
            encoder -> leftHaps[i] = (LEFT_ALLELE(encoder -> genos[i]) + 1) * LOG_PRIMES[0];
            encoder -> rightHaps[i] = (RIGHT_ALLELE(encoder -> genos[i]) + 1) * LOG_PRIMES[0];

            if (collapseMissingGenotypes && (LEFT_ALLELE(encoder -> genos[i]) == numAlleles || RIGHT_ALLELE(encoder -> genos[i]) == numAlleles)) {
                encoder -> leftHaps[i] = MISSING;
                encoder -> rightHaps[i] = MISSING;
            }
        } else if (collapseMissingGenotypes && (encoder -> leftHaps[i] == MISSING || encoder -> rightHaps[i] == MISSING || LEFT_ALLELE(encoder -> genos[i]) == numAlleles || RIGHT_ALLELE(encoder -> genos[i]) == numAlleles)) {
            encoder -> leftHaps[i] = MISSING;
            encoder -> rightHaps[i] = MISSING;
        } else {
            encoder -> leftHaps[i] += (LEFT_ALLELE(encoder -> genos[i]) + 1) * LOG_PRIMES[encoder -> numLoci];
            encoder -> rightHaps[i] += (RIGHT_ALLELE(encoder -> genos[i]) + 1) * LOG_PRIMES[encoder -> numLoci];
        }
    }

    encoder -> numLoci++;

}

bool get_next_haplotype(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, bool collapseMissingGenotypes, int HAP_SIZE) {

    // If EOF, there is no next haplotype.
    if (parser -> isEOF)
        return false;
    
    // Set chromosome of next record from parser.
    encoder -> chrom -> l = 0;
    kputs(ks_str(parser -> nextChrom), encoder -> chrom);

    // Set start locus of next record from parser.
    encoder -> startLocus = parser -> nextPos;

    // Empty haplotype.
    encoder -> numLoci = 0;

    // Used to flag if haplotype is on the same chromosome.
    bool isSameChrom = true;

    // Holds number of alleles at each record.
    int numAlleles;

    // Create the haplotype.
    while(!(parser -> isEOF) && (encoder -> numLoci < HAP_SIZE) && isSameChrom) {
        // Get the next record from the VCF file.
        get_next_locus(parser, encoder -> chrom, &(encoder -> endLocus), &numAlleles, &(encoder -> genos));
        // Add locus to haplotype.
        add_locus(encoder, numAlleles, collapseMissingGenotypes);
        // Make sure the next locus is on the same haplotype.
        isSameChrom = strcmp(ks_str(encoder -> chrom), ks_str(parser -> nextChrom)) == 0;
    }

    // Not EOF, complete haplotype, and next loci is on the same chromsome.
    return !(parser -> isEOF) && encoder -> numLoci == HAP_SIZE && isSameChrom;

}

void destroy_haplotype_encoder(HaplotypeEncoder* encoder) {

    // Free memory used by arrays.
    free(encoder -> genos);
    free(encoder -> leftHaps);
    free(encoder -> rightHaps);
    // Free chromosome string.
    free(ks_str(encoder -> chrom)); free(encoder -> chrom);
    // Free structure.
    free(encoder);

}

double LOG_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541
};


// Used to test the encoder.

/*
void print_encoder_info(HaplotypeEncoder* encoder) {
    printf("Chromosome: %s\n", ks_str(encoder -> chrom));
    printf("Start locus: %d\n", encoder -> startLocus);
    printf("End locus: %d\n", encoder -> endLocus);
    printf("Number of loci: %d\n", encoder -> numLoci);
    printf("Sample Haplotypes:\n");
    for (int i = 0; i < encoder -> numSamples; i++)
        printf("Sample %d -> %lf/%lf\n", i + 1, encoder -> leftHaps[i], encoder -> rightHaps[i]);
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
*/