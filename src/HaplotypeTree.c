
#include "HaplotypeTree.h"

#define LEFT_ALLELE(a) (a >> 4)
#define RIGHT_ALLELE(a) (a & 0x0F)

HaplotypeTree* init_haplotype_tree(int numSamples) {

    HaplotypeTree* tree = (HaplotypeTree*) calloc(1, sizeof(HaplotypeTree));

    tree -> numSamples = numSamples;
    tree -> genotypes = (GENOTYPE*) calloc(numSamples, sizeof(GENOTYPE));
    tree -> leftHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));
    tree -> rightHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));

    tree -> chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));

    tree -> labelMap = kh_init(32);

    tree -> numLeaves = 1;

    return tree;

}

void add_locus(HaplotypeTree* tree, int numAlleles, bool collapseMissingGenotypes) {
    
    for (int i = 0; i < tree -> numSamples; i++) {
        if (tree -> numLeaves == 1) {
            tree -> leftHaplotype[i] = LEFT_ALLELE(tree -> genotypes[i]);
            tree -> rightHaplotype[i] = RIGHT_ALLELE(tree -> genotypes[i]);
            if (collapseMissingGenotypes && (tree -> leftHaplotype[i] == numAlleles || tree -> rightHaplotype[i] == numAlleles)) {
                tree -> leftHaplotype[i] = numAlleles;
                tree -> rightHaplotype[i] = numAlleles;
            }
        } else if (collapseMissingGenotypes && (tree -> leftHaplotype[i] == (tree -> numLeaves - 1) || tree -> rightHaplotype[i] == (tree -> numLeaves - 1) || LEFT_ALLELE(tree -> genotypes[i]) == numAlleles || RIGHT_ALLELE(tree -> genotypes[i]) == numAlleles)) {
            tree -> leftHaplotype[i] = tree -> numLeaves * (numAlleles + 1) - 1;
            tree -> rightHaplotype[i] = tree -> numLeaves * (numAlleles + 1) - 1;
        } else {
            tree -> leftHaplotype[i] = tree -> leftHaplotype[i] * (numAlleles + 1) + LEFT_ALLELE(tree -> genotypes[i]);
            tree -> rightHaplotype[i] = tree -> rightHaplotype[i] * (numAlleles + 1) + RIGHT_ALLELE(tree -> genotypes[i]);
        }
    }
    
    tree -> numLeaves = (tree -> numLeaves) * (numAlleles + 1);

    if (tree -> numLeaves >= MAX_NUM_LEAVES)
        relabel_haplotypes(tree);

}

void relabel_haplotypes(HaplotypeTree* tree) {

    kh_clear(32, tree -> labelMap);
    
    int ret, newLabel = 0;

    khiter_t k = kh_put(32, tree -> labelMap, tree -> numLeaves - 1, &ret);
    kh_value(tree -> labelMap, k) = 0xFFFFFFFF;

    for (int i = 0; i < tree -> numSamples; i++) {

        if (kh_get(32, tree -> labelMap, tree -> leftHaplotype[i]) == kh_end(tree -> labelMap)) {
            k = kh_put(32, tree -> labelMap, tree -> leftHaplotype[i], &ret);
            kh_value(tree -> labelMap, k) = newLabel++;
        }

        if (kh_get(32, tree -> labelMap, tree -> rightHaplotype[i]) == kh_end(tree -> labelMap)) {
            k = kh_put(32, tree -> labelMap, tree -> rightHaplotype[i], &ret);
            kh_value(tree -> labelMap, k) = newLabel++;
        }
        
        tree -> leftHaplotype[i] = kh_value(tree -> labelMap, kh_get(32, tree -> labelMap, tree -> leftHaplotype[i]));
        tree -> rightHaplotype[i] = kh_value(tree -> labelMap, kh_get(32, tree -> labelMap, tree -> rightHaplotype[i]));

    }

    for (int i = 0; i < tree -> numSamples; i++) {
        if (tree -> leftHaplotype[i] == 0xFFFFFFFF)
            tree -> leftHaplotype[i] = newLabel;
        if (tree -> rightHaplotype[i] == 0xFFFFFFFF)
            tree -> rightHaplotype[i] = newLabel;
    }

    tree -> numLeaves = newLabel + 1;

}

bool get_next_haplotype(VCFGenotypeParser* parser, HaplotypeTree* tree, int HAP_SIZE) {

    if (parser -> isEOF)
        return false;

    tree -> chromosome -> l = 0;
    kputs(ks_str(parser -> nextChromosome), tree -> chromosome);
    tree -> startLocus = parser -> nextPosition;

    tree -> numLeaves = 1;
    tree -> numLoci = 0;

    bool isSameChromosome = true;
    int numAlleles;

    while(!(parser -> isEOF) && (tree -> numLoci < HAP_SIZE) && isSameChromosome) {
        get_next_locus(parser, tree -> chromosome, &(tree -> endLocus), &numAlleles, &(tree -> genotypes));
        add_locus(tree, numAlleles, true);
        isSameChromosome = strcmp(ks_str(tree -> chromosome), ks_str(parser -> nextChromosome)) == 0;
        tree -> numLoci++;
    }

    return tree -> numLoci == HAP_SIZE && isSameChromosome;

}

void destroy_haplotype_tree(HaplotypeTree* tree) {

    free(tree -> genotypes);
    free(tree -> leftHaplotype);
    free(tree -> rightHaplotype);
    free(ks_str(tree -> chromosome)); free(tree -> chromosome);
    kh_destroy(32, tree -> labelMap);
    free(tree);

}

/*
void print_tree_info(HaplotypeTree* tree) {
    printf("Chromosome: %s\n", ks_str(tree -> chromosome));
    printf("Start locus: %d\n", tree -> startLocus);
    printf("End locus: %d\n", tree -> endLocus);
    printf("Number of loci: %d\n", tree -> numLoci);
    printf("Sample Haplotypes:\n");
    for (int i = 0; i < tree -> numSamples; i++)
        printf("Sample %d -> %d/%d\n", i + 1, tree -> leftHaplotype[i], tree -> rightHaplotype[i]);
}

int main() {

    VCFGenotypeParser* parser = init_vcf_genotype_parser("haplotype_tree_test.vcf.gz");
    HaplotypeTree* tree = init_haplotype_tree(parser -> num_samples);
    printf("\nTest 1\n");
    printf("---------\n");
    printf("Read in whole chromosomes:\n\n");
    get_next_haplotype(parser, tree, 10);
    printf("First Haplotype:\n");
    print_tree_info(tree);
    get_next_haplotype(parser, tree, 10);
    printf("\nSecond Haplotype:\n");
    print_tree_info(tree);
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    parser = init_vcf_genotype_parser("haplotype_tree_test.vcf.gz");
    printf("\nTest 2\n");
    printf("---------\n");
    printf("Read in 2-loci haplotypes:\n\n");
    get_next_haplotype(parser, tree, 2);
    printf("First Haplotype:\n");
    print_tree_info(tree);
    get_next_haplotype(parser, tree, 2);
    printf("\nSecond Haplotype:\n");
    print_tree_info(tree);
    get_next_haplotype(parser, tree, 2);
    printf("\nThird Haplotype:\n");
    print_tree_info(tree);
    destroy_vcf_genotype_parser(parser);

    printf("\n");
    parser = init_vcf_genotype_parser("haplotype_tree_test.vcf.gz");
    printf("\nTest 3\n");
    printf("---------\n");
    printf("Read in 1-locus haplotypes:\n\n");
    get_next_haplotype(parser, tree, 1);
    printf("First Haplotype:\n");
    print_tree_info(tree);
    get_next_haplotype(parser, tree, 1);
    printf("\nSecond Haplotype:\n");
    print_tree_info(tree);
    get_next_haplotype(parser, tree, 1);
    printf("\nThird Haplotype:\n");
    print_tree_info(tree);
    get_next_haplotype(parser, tree, 1);
    printf("\nFourth Haplotype:\n");
    print_tree_info(tree);
    get_next_haplotype(parser, tree, 1);
    printf("\nFifth Haplotype:\n");
    print_tree_info(tree);
    destroy_vcf_genotype_parser(parser);

    printf("\n");

    destroy_haplotype_tree(tree);

}
*/