
// File: HaplotypeEncoder.c
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Track the haplotypes of each sample using a simplified arithmetic encoding.

#include "HaplotypeEncoder.h"

// Macros to get the left and right allele from the 1-byte encoding.
#define LEFT_ALLELE(a) (a >> 4)
#define RIGHT_ALLELE(a) (a & 0x0F)

HaplotypeEncoder* init_haplotype_encoder(int numSamples) {

    // Allocate the structure's memory.
    HaplotypeEncoder* encoder = (HaplotypeEncoder*) calloc(1, sizeof(HaplotypeEncoder));

    // Allocate memory used for arrays.
    encoder -> numSamples = numSamples;
    encoder -> genotypes = (GENOTYPE*) calloc(numSamples, sizeof(GENOTYPE));
    encoder -> leftHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));
    encoder -> rightHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));

    // Allocate string to hold chromosome.
    encoder -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));

    // Create the hash table used for relabeling.
    encoder -> labelMap = kh_init(32);
    
    // The tree starts off with one leaf, the empty string.
    encoder -> numLeaves = 1;

    // Return the encoder.
    return encoder;

}

void add_locus(HaplotypeEncoder* encoder, int numAlleles, bool collapseMissingGenotypes) {
    
    // Iterate through the samples.
    for (int i = 0; i < encoder -> numSamples; i++) {
        // Haplotypes get the allele genotypes at the first level.
        if (encoder -> numLeaves == 1) {
            encoder -> leftHaplotype[i] = LEFT_ALLELE(encoder -> genotypes[i]);
            encoder -> rightHaplotype[i] = RIGHT_ALLELE(encoder -> genotypes[i]);
            // If we are collapsing missing genotypes, and either allele is missing, move both to the right most leaf.
            if (collapseMissingGenotypes && (encoder -> leftHaplotype[i] == numAlleles || encoder -> rightHaplotype[i] == numAlleles)) {
                encoder -> leftHaplotype[i] = numAlleles;
                encoder -> rightHaplotype[i] = numAlleles;
            }
        // If we are collapsing genotypes and a missing genotype is encountered, move each to the left most leaf.
        } else if (collapseMissingGenotypes && (encoder -> leftHaplotype[i] == (encoder -> numLeaves - 1) || encoder -> rightHaplotype[i] == (encoder -> numLeaves - 1) || LEFT_ALLELE(encoder -> genotypes[i]) == numAlleles || RIGHT_ALLELE(encoder -> genotypes[i]) == numAlleles)) {
            encoder -> leftHaplotype[i] = encoder -> numLeaves * (numAlleles + 1) - 1;
            encoder -> rightHaplotype[i] = encoder -> numLeaves * (numAlleles + 1) - 1;
        // Otherwise, advance haplotypes to the next level.
        } else {
            encoder -> leftHaplotype[i] = encoder -> leftHaplotype[i] * (numAlleles + 1) + LEFT_ALLELE(encoder -> genotypes[i]);
            encoder -> rightHaplotype[i] = encoder -> rightHaplotype[i] * (numAlleles + 1) + RIGHT_ALLELE(encoder -> genotypes[i]);
        }
    }
    
    // Extend tree.
    encoder -> numLeaves = (encoder -> numLeaves) * (numAlleles + 1);

    // If max number of leaves is succeeded, then relabel tree.
    //  This will create a tree with a maximum of 2 * numSamples leaves.
    if (encoder -> numLeaves >= MAX_NUM_LEAVES)
        relabel_haplotypes(encoder);

}

void relabel_haplotypes(HaplotypeEncoder* encoder) {

    // Clear hash table of any contents without deallocating memory.
    kh_clear(32, encoder -> labelMap);
    
    // Track the new label.
    int ret, newLabel = 0;

    // Map the right most leaf to 0xFFFFFFFF.
    khiter_t k = kh_put(32, encoder -> labelMap, encoder -> numLeaves - 1, &ret);
    kh_value(encoder -> labelMap, k) = 0xFFFFFFFF;

    for (int i = 0; i < encoder -> numSamples; i++) {

        // If encoded value does not exist in the hash table, map value to new label.
        if (kh_get(32, encoder -> labelMap, encoder -> leftHaplotype[i]) == kh_end(encoder -> labelMap)) {
            k = kh_put(32, encoder -> labelMap, encoder -> leftHaplotype[i], &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }

        // If encoded value does not exist in the hash table, map value to new label.
        if (kh_get(32, encoder -> labelMap, encoder -> rightHaplotype[i]) == kh_end(encoder -> labelMap)) {
            k = kh_put(32, encoder -> labelMap, encoder -> rightHaplotype[i], &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        
        // Relabel haplotypes.
        encoder -> leftHaplotype[i] = kh_value(encoder -> labelMap, kh_get(32, encoder -> labelMap, encoder -> leftHaplotype[i]));
        encoder -> rightHaplotype[i] = kh_value(encoder -> labelMap, kh_get(32, encoder -> labelMap, encoder -> rightHaplotype[i]));

    }

    // For all of the haplotypes labeled with the right-most leaf of old tree, assign right-most label of new tree.
    for (int i = 0; i < encoder -> numSamples; i++) {
        if (encoder -> leftHaplotype[i] == 0xFFFFFFFF)
            encoder -> leftHaplotype[i] = newLabel;
        if (encoder -> rightHaplotype[i] == 0xFFFFFFFF)
            encoder -> rightHaplotype[i] = newLabel;
    }

    // New number of leaves.
    encoder -> numLeaves = newLabel + 1;

}

bool get_next_haplotype(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, bool collapseMissingGenotypes, int HAP_SIZE) {

    // If EOF, there is no next haplotype.
    if (parser -> isEOF)
        return false;
    
    // Set chromosome of next record from parser.
    encoder -> chromosome -> l = 0;
    kputs(ks_str(parser -> nextChromosome), encoder -> chromosome);

    // Set start locus of next record from parser.
    encoder -> startLocus = parser -> nextPosition;

    // Reset tree.
    encoder -> numLeaves = 1;

    // Empty haplotype.
    encoder -> numLoci = 0;

    // Used to flag if haplotype is on the same chromosome.
    bool isSameChromosome = true;

    // Holds number of alleles at each record.
    int numAlleles;

    // Create the haplotype.
    while(!(parser -> isEOF) && (encoder -> numLoci < HAP_SIZE) && isSameChromosome) {
        // Get the next record from the VCF file.
        get_next_locus(parser, encoder -> chromosome, &(encoder -> endLocus), &numAlleles, &(encoder -> genotypes));
        // Add locus to haplotype.
        add_locus(encoder, numAlleles, collapseMissingGenotypes);
        // Make sure the next locus is on the same haplotype.
        isSameChromosome = strcmp(ks_str(encoder -> chromosome), ks_str(parser -> nextChromosome)) == 0;
        encoder -> numLoci++;
    }

    // Not EOF, complete haplotype, and next loci is on the same chromsome.
    return !(parser -> isEOF) && encoder -> numLoci == HAP_SIZE && isSameChromosome;

}

void destroy_haplotype_encoder(HaplotypeEncoder* encoder) {

    // Free memory used by arrays.
    free(encoder -> genotypes);
    free(encoder -> leftHaplotype);
    free(encoder -> rightHaplotype);
    // Free chromosome string.
    free(ks_str(encoder -> chromosome)); free(encoder -> chromosome);
    // Free hash map.
    kh_destroy(32, encoder -> labelMap);
    // Free structure.
    free(encoder);

}

// Used to test the haplotype encoder.

/*
void print_encoder_info(HaplotypeEncoder* encoder) {
    printf("Chromosome: %s\n", ks_str(encoder -> chromosome));
    printf("Start locus: %d\n", encoder -> startLocus);
    printf("End locus: %d\n", encoder -> endLocus);
    printf("Number of loci: %d\n", encoder -> numLoci);
    printf("Sample Haplotypes:\n");
    for (int i = 0; i < encoder -> numSamples; i++)
        printf("Sample %d -> %d/%d\n", i + 1, encoder -> leftHaplotype[i], encoder -> rightHaplotype[i]);
}

int main() {

    VCFGenotypeParser* parser = init_vcf_genotype_parser("haplotype_encoder_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> num_samples);
    printf("\nTest 1\n");
    printf("---------\n");
    printf("Read in whole chromosomes:\n\n");
    get_next_haplotype(parser, encoder, 10);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 10);
    printf("\nSecond Haplotype:\n");
    print_encoder_info(encoder);
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    parser = init_vcf_genotype_parser("haplotype_encoder_test.vcf.gz");
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
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    parser = init_vcf_genotype_parser("haplotype_encoder_test.vcf.gz");
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
    destroy_vcf_genotype_parser(parser);

    printf("\n");

    destroy_haplotype_encoder(encoder);

}
*/