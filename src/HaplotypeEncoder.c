
// File: HaplotypeEncoder.c
// Date: 6 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Assume a VCF file is phased and encode haplotypes of samples' genotypes.

#include "HaplotypeEncoder.h"

HaplotypeEncoder_t* init_haplotype_encoder(int numSamples) {

    // Create structure.
    HaplotypeEncoder_t* encoder = (HaplotypeEncoder_t*) calloc(1, sizeof(HaplotypeEncoder_t));

    // Allocate reguired memory.
    encoder -> numSamples = numSamples;
    encoder -> locus = (Locus*) calloc(numSamples, sizeof(Locus));
    encoder -> genotypes = (Genotype_t*) calloc(numSamples, sizeof(Genotype_t));
    encoder -> chrom = init_kstring(NULL);
    encoder -> labelMap = kh_init(haplotype);
    // The tree has one node, which corresponds to the empty string.
    encoder -> numLeaves = 1;

    return encoder;

}

// The method used to relabel haplotypes when the tree gets too large.
// Accepts:
//  HaplotypeEncoder_t* encoder -> The encoder whose haplotypes we are going to relabel.
// Returns: void.
void relabel_haplotypes(HaplotypeEncoder_t* encoder) {

    khint_t j;
    khiter_t k;
    int ret;

    // Before relabeling, we mark all the current labels in the hash table with a flag to signal that
    //  they can be reassigned. I debated if this approach is better than clearing buckets and repopulating
    //  the hash table. This maybe faster than reallocating memory repeatedly.
    for (j = kh_begin(encoder -> labelMap); j != kh_end(encoder -> labelMap); j++) {
		if (!kh_exist(encoder -> labelMap, j)) 
            continue;
		kh_val(encoder -> labelMap, j) = MISSING;
	}
    
    // We start the new label at 0.
    Haplotype newLabel = 0;

    // For each of the samples ...
    for (int i = 0; i < encoder -> numSamples; i++) {
        // If the sample has missing genotypes, we skip relabeling.
        if (encoder -> genotypes[i].left == MISSING)
            continue;
        
        // Label left haplotype.
        k = kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left);
        // If the left haplotype encoding has never been encountered before ...
        if (k == kh_end(encoder -> labelMap)) {
            // Insert encoding and generate label for the haplotype.
            k = kh_put(haplotype, encoder -> labelMap, encoder -> genotypes[i].left, &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        // If the encoding has been used before in a previous relabeling but needs a new label ...
        if (kh_value(encoder -> labelMap, k) == MISSING) {
            // Create new label.
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        // Get label from hash table and relabel left haplotype.
        encoder -> genotypes[i].left = kh_value(encoder -> labelMap, kh_get(haplotype, encoder -> labelMap, encoder -> genotypes[i].left));

        // Do the same for the right haplotype as the left.
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

    // The number of leaves in the tree is the same as the number of new labels (unique haplotypes).
    encoder -> numLeaves = newLabel;

}

// Add a new locus to each haplotype and update encodings.
//  We are extending our tree.
// Accepts:
//  HaplotypeEncoder_t* encoder -> The encoder with the haplotypes we are extending.
//  int numAlleles -> The number of alleles at the locus we are appending to the haplotypes.
// Returns: void.
void add_locus(HaplotypeEncoder_t* encoder, int numAlleles) {

    // If extending our tree causes integer overflow, we relabel the haplotypes.
    if (numAlleles * encoder -> numLeaves == MISSING || (numAlleles * encoder -> numLeaves) / numAlleles != encoder -> numLeaves)
        relabel_haplotypes(encoder);

    // NOTE: We are assuming that the new locus is in encoder -> locus.
    // For each of the samples ...
    for (int i = 0; i < encoder -> numSamples; i++) {
        // If this locus is the first in the tree.
        if (encoder -> numLeaves == 1) {
            // Assign the left and right alleles to their respective haplotypes.
            encoder -> genotypes[i].left = LEFT_ALLELE(encoder -> locus[i]);
            encoder -> genotypes[i].right = RIGHT_ALLELE(encoder -> locus[i]);
            // If either of the alleles are the missing genotype, both the left and right are flagged as missing.
            if (encoder -> genotypes[i].left == numAlleles || encoder -> genotypes[i].right == numAlleles) {
                encoder -> genotypes[i].left = MISSING;
                encoder -> genotypes[i].right = MISSING;
            }
        // If either the left and right haplotypes are missing or the current loci contains a missing genotype, we flag both haplotypes as missing. 
        } else if (encoder -> genotypes[i].left == MISSING || LEFT_ALLELE(encoder -> locus[i]) == numAlleles || RIGHT_ALLELE(encoder -> locus[i]) == numAlleles) {
            encoder -> genotypes[i].left = MISSING;
            encoder -> genotypes[i].right = MISSING;
        // Otherwise, both genotypes are not missing, so we move the haplotypes to their new leaves.
        } else {
            encoder -> genotypes[i].left = encoder -> genotypes[i].left * numAlleles + LEFT_ALLELE(encoder -> locus[i]);
            encoder -> genotypes[i].right = encoder -> genotypes[i].right * numAlleles + RIGHT_ALLELE(encoder -> locus[i]);
        }
    }
    
    // We extend the tree.
    encoder -> numLeaves = (encoder -> numLeaves) * numAlleles;

}

bool get_next_haplotype(VCFLocusParser_t* parser, HaplotypeEncoder_t* encoder, int HAP_SIZE) {

    // If the end of the VCF file has been reached, we cannot get another haplotype.
    if (parser -> isEOF)
        return false;
    
    // Set chromosome and start coordinate of the haplotype.
    ks_overwrite(ks_str(parser -> nextChrom), encoder -> chrom);
    encoder -> startCoord = parser -> nextCoord;

    // Reset the tree.
    encoder -> numLeaves = 1;

    // There are no loci in the haplotype yet.
    encoder -> numLoci = 0;

    // Used to test if the next locus is on the same coordinate as the current.
    bool isSameChrom = true;

    // Holds the number of alleles at a locus.
    int numAlleles;

    // While the end of VCF file has not been reached the maximum haplotype size has not been reached
    //  and the current locus is on the same chromosome as the next locus.
    while(!(parser -> isEOF) && (encoder -> numLoci < HAP_SIZE) && isSameChrom) {
        // Get the next locus from the VCF file.
        get_next_locus(parser, encoder -> chrom, &(encoder -> endCoord), &numAlleles, &(encoder -> locus));
        // Extend the haplotypes.
        add_locus(encoder, numAlleles);
        isSameChrom = strcmp(ks_str(encoder -> chrom), ks_str(parser -> nextChrom)) == 0;
        encoder -> numLoci++;
    }

    return !(parser -> isEOF) && encoder -> numLoci == HAP_SIZE && isSameChrom;

}

void destroy_haplotype_encoder(HaplotypeEncoder_t* encoder) {
    if (encoder == NULL)
        return;
    free(encoder -> locus);
    free(encoder -> genotypes);
    destroy_kstring(encoder -> chrom);
    kh_destroy(haplotype, encoder -> labelMap);
    free(encoder);
}

// Used to test the haplotype encoder.

void print_encoder_info(HaplotypeEncoder_t* encoder) {
    printf("Chromosome: %s\n", ks_str(encoder -> chrom));
    printf("Start locus: %d\n", encoder -> startCoord);
    printf("End locus: %d\n", encoder -> endCoord);
    printf("Number of loci: %d\n", encoder -> numLoci);
    printf("Sample Haplotypes:\n");
    for (int i = 0; i < encoder -> numSamples; i++)
        printf("Sample %d -> %ld/%ld\n", i + 1, encoder -> genotypes[i].left, encoder -> genotypes[i].right);
}

/*
int main() {

    VCFLocusParser_t* parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz", NULL, false, 0, 1, false);
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(parser -> numSamples);
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
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz", NULL, false, 0, 1, false);
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
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test.vcf.gz", NULL, false, 0, 1, false);
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
    parser = init_vcf_locus_parser("./data/haplotype_encoder_test2.vcf.gz", NULL, false, 0, 1, false);
    printf("\nTest 4\n");
    printf("---------\n");
    get_next_haplotype(parser, encoder, 4);
    printf("First Haplotype:\n");
    print_encoder_info(encoder);
    printf("\nRelabeled:\n");
    relabel_haplotypes(encoder);
    print_encoder_info(encoder);
    get_next_haplotype(parser, encoder, 10);
    printf("\nRe-labeled the Relabels:\n");
    print_encoder_info(encoder);
    printf("\nRelabeled:\n");
    relabel_haplotypes(encoder);
    print_encoder_info(encoder);
    destroy_vcf_locus_parser(parser);

    printf("\n");
    
    destroy_haplotype_encoder(encoder);
}
*/