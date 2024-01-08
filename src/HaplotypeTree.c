
#include "HaplotypeTree.h"

#define LEFT_ALLELE(a) (a >> 4)
#define RIGHT_ALLELE(a) (a & 0x0F)

HaplotypeTree* init_haplotype_tree(int numSamples) {

    HaplotypeTree* tree = (HaplotypeTree*) calloc(1, sizeof(HaplotypeTree));

    tree -> numSamples = numSamples;
    tree -> leftHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));
    tree -> rightHaplotype = (unsigned int*) calloc(numSamples, sizeof(int));

    tree -> labelMap = kh_init(32);

    tree -> numLeaves = 1;

    return tree;

}

void add_locus(HaplotypeTree* tree, int numAlleles, GENOTYPE* genotypes, bool collapseMissingGenotypes) {
    
    for (int i = 0; i < tree -> numSamples; i++) {
        if (tree -> numLeaves == 1) {
            tree -> leftHaplotype[i] = LEFT_ALLELE(genotypes[i]);
            tree -> rightHaplotype[i] = RIGHT_ALLELE(genotypes[i]);
        } else if (collapseMissingGenotypes && (tree -> leftHaplotype[i] == (tree -> numLeaves - 1) || tree -> rightHaplotype[i] == (tree -> numLeaves - 1) || LEFT_ALLELE(genotypes[i]) == numAlleles || RIGHT_ALLELE(genotypes[i]) == numAlleles)) {
            tree -> leftHaplotype[i] = tree -> numLeaves * (numAlleles + 1) - 1;
            tree -> rightHaplotype[i] = tree -> numLeaves * (numAlleles + 1) - 1;
        } else {
            tree -> leftHaplotype[i] = tree -> leftHaplotype[i] * (numAlleles + 1) + LEFT_ALLELE(genotypes[i]);
            tree -> rightHaplotype[i] = tree -> rightHaplotype[i] * (numAlleles + 1) + RIGHT_ALLELE(genotypes[i]);
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

void destroy_haplotype_tree(HaplotypeTree* tree) {

    free(tree -> leftHaplotype);
    free(tree -> rightHaplotype);
    kh_destroy(32, tree -> labelMap);
    free(tree);

}

int main() {

    VCFGenotypeParser* parser = init_vcf_genotype_parser("haplotype_tree_test.vcf.gz");
    
    kstring_t* chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    int position;
    int numOfAlleles;
    GENOTYPE* genotypes = (GENOTYPE*) calloc(parser -> num_samples, sizeof(GENOTYPE));

    HaplotypeTree* tree = init_haplotype_tree(parser -> num_samples);

    

    destroy_vcf_genotype_parser(parser);
    free(genotypes);
    free(ks_str(chromosome)); free(chromosome);

}