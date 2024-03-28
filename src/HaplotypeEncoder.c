
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
    encoder -> labelMap = kh_init(64);

    return encoder;

}

void relabel_haplotypes(HaplotypeEncoder* encoder) {

    kh_clear(64, encoder -> labelMap);
    
    int ret, newLabel = 0;

    khiter_t k = kh_put(64, encoder -> labelMap, MISSING, &ret);
    kh_value(encoder -> labelMap, k) = MISSING;

    Haplotype left, right;

    for (int i = 0; i < encoder -> numSamples; i++) {

        left = encoder -> genotypes[i].left;
        right = encoder -> genotypes[i].right;

        if (kh_get(64, encoder -> labelMap, left) == kh_end(encoder -> labelMap)) {
            k = kh_put(64, encoder -> labelMap, left, &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }

        if (kh_get(64, encoder -> labelMap, right) == kh_end(encoder -> labelMap)) {
            k = kh_put(64, encoder -> labelMap, right, &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        
        encoder -> genotypes[i].left = kh_value(encoder -> labelMap, kh_get(64, encoder -> labelMap, left));
        encoder -> genotypes[i].right = kh_value(encoder -> labelMap, kh_get(64, encoder -> labelMap, right));
    }

    encoder -> numLeaves = newLabel + 1;

}

void add_locus(HaplotypeEncoder* encoder, int numAlleles, bool collapseMissingGenotypes) {

    if (encoder -> numLeaves * numAlleles == MISSING || (encoder -> numLeaves * numAlleles) / numAlleles != encoder -> numLeaves)
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
            encoder -> genotypes[i].left = encoder -> genotypes[i].left * (numAlleles + 1) + LEFT_ALLELE(encoder -> locus[i]);
            encoder -> genotypes[i].right = encoder -> genotypes[i].right * (numAlleles + 1) + RIGHT_ALLELE(encoder -> locus[i]);
        }
    }
    
    encoder -> numLeaves = (encoder -> numLeaves) * (numAlleles + 1);

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
    }

    return !(parser -> isEOF) && encoder -> numLoci == HAP_SIZE && isSameChrom;

}

void destroy_haplotype_encoder(HaplotypeEncoder* encoder) {
    if (encoder == NULL)
        return;
    free(encoder -> locus);
    free(encoder -> genotypes);
    free(ks_str(encoder -> chrom)); free(encoder -> chrom);
    kh_destroy(64, encoder -> labelMap);
    free(encoder);
}