
#include "HaplotypeEncoder.h"

#define LEFT_ALLELE(a) (a >> 4)
#define RIGHT_ALLELE(a) (a & 0x0F)

HaplotypeEncoder* init_haplotype_encoder(int numSamples) {

    HaplotypeEncoder* encoder = (HaplotypeEncoder*) calloc(1, sizeof(HaplotypeEncoder));

    encoder -> numSamples = numSamples;
    encoder -> genotypes = (GENOTYPE*) calloc(numSamples, sizeof(GENOTYPE));
    encoder -> leftHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));
    encoder -> rightHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));

    encoder -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));

    encoder -> labelMap = kh_init(32);

    encoder -> numLeaves = 1;

    return encoder;

}

void add_locus(HaplotypeEncoder* encoder, int numAlleles, bool collapseMissingGenotypes) {
    
    for (int i = 0; i < encoder -> numSamples; i++) {
        if (encoder -> numLeaves == 1) {
            encoder -> leftHaplotype[i] = LEFT_ALLELE(encoder -> genotypes[i]);
            encoder -> rightHaplotype[i] = RIGHT_ALLELE(encoder -> genotypes[i]);
            if (collapseMissingGenotypes && (encoder -> leftHaplotype[i] == numAlleles || encoder -> rightHaplotype[i] == numAlleles)) {
                encoder -> leftHaplotype[i] = numAlleles;
                encoder -> rightHaplotype[i] = numAlleles;
            }
        } else if (collapseMissingGenotypes && (encoder -> leftHaplotype[i] == (encoder -> numLeaves - 1) || encoder -> rightHaplotype[i] == (encoder -> numLeaves - 1) || LEFT_ALLELE(encoder -> genotypes[i]) == numAlleles || RIGHT_ALLELE(encoder -> genotypes[i]) == numAlleles)) {
            encoder -> leftHaplotype[i] = encoder -> numLeaves * (numAlleles + 1) - 1;
            encoder -> rightHaplotype[i] = encoder -> numLeaves * (numAlleles + 1) - 1;
        } else {
            encoder -> leftHaplotype[i] = encoder -> leftHaplotype[i] * (numAlleles + 1) + LEFT_ALLELE(encoder -> genotypes[i]);
            encoder -> rightHaplotype[i] = encoder -> rightHaplotype[i] * (numAlleles + 1) + RIGHT_ALLELE(encoder -> genotypes[i]);
        }
    }
    
    encoder -> numLeaves = (encoder -> numLeaves) * (numAlleles + 1);

    if (encoder -> numLeaves >= MAX_NUM_LEAVES)
        relabel_haplotypes(encoder);

}

void relabel_haplotypes(HaplotypeEncoder* encoder) {

    kh_clear(32, encoder -> labelMap);
    
    int ret, newLabel = 0;

    khiter_t k = kh_put(32, encoder -> labelMap, encoder -> numLeaves - 1, &ret);
    kh_value(encoder -> labelMap, k) = 0xFFFFFFFF;

    for (int i = 0; i < encoder -> numSamples; i++) {

        if (kh_get(32, encoder -> labelMap, encoder -> leftHaplotype[i]) == kh_end(encoder -> labelMap)) {
            k = kh_put(32, encoder -> labelMap, encoder -> leftHaplotype[i], &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }

        if (kh_get(32, encoder -> labelMap, encoder -> rightHaplotype[i]) == kh_end(encoder -> labelMap)) {
            k = kh_put(32, encoder -> labelMap, encoder -> rightHaplotype[i], &ret);
            kh_value(encoder -> labelMap, k) = newLabel++;
        }
        
        encoder -> leftHaplotype[i] = kh_value(encoder -> labelMap, kh_get(32, encoder -> labelMap, encoder -> leftHaplotype[i]));
        encoder -> rightHaplotype[i] = kh_value(encoder -> labelMap, kh_get(32, encoder -> labelMap, encoder -> rightHaplotype[i]));

    }

    for (int i = 0; i < encoder -> numSamples; i++) {
        if (encoder -> leftHaplotype[i] == 0xFFFFFFFF)
            encoder -> leftHaplotype[i] = newLabel;
        if (encoder -> rightHaplotype[i] == 0xFFFFFFFF)
            encoder -> rightHaplotype[i] = newLabel;
    }

    encoder -> numLeaves = newLabel + 1;

}

bool get_next_haplotype(VCFGenotypeParser* parser, HaplotypeEncoder* encoder, int HAP_SIZE) {

    if (parser -> isEOF)
        return false;

    encoder -> chromosome -> l = 0;
    kputs(ks_str(parser -> nextChromosome), encoder -> chromosome);
    encoder -> startLocus = parser -> nextPosition;

    encoder -> numLeaves = 1;
    encoder -> numLoci = 0;

    bool isSameChromosome = true;
    int numAlleles;

    while(!(parser -> isEOF) && (encoder -> numLoci < HAP_SIZE) && isSameChromosome) {
        get_next_locus(parser, encoder -> chromosome, &(encoder -> endLocus), &numAlleles, &(encoder -> genotypes));
        add_locus(encoder, numAlleles, true);
        isSameChromosome = strcmp(ks_str(encoder -> chromosome), ks_str(parser -> nextChromosome)) == 0;
        encoder -> numLoci++;
    }

    return encoder -> numLoci == HAP_SIZE && isSameChromosome;

}

void destroy_haplotype_encoder(HaplotypeEncoder* encoder) {

    free(encoder -> genotypes);
    free(encoder -> leftHaplotype);
    free(encoder -> rightHaplotype);
    free(ks_str(encoder -> chromosome)); free(encoder -> chromosome);
    kh_destroy(32, encoder -> labelMap);
    free(encoder);

}

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