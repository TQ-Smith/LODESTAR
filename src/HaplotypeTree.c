
#include "HaplotypeTree.h"

#define LEFT_ALLELE(a) (a >> 4)
#define RIGHT_ALLELE(a) (a & 0x0F)

#define CALCULATE_NEW_LABEL(curLabel, numNodes, prevGroups, prevAlleles, numAlleles, allele) \
    (curLabel + (numNodes - (prevGroups * prevAlleles) - 1) * numAlleles + (numNodes - curLabel - 1) + allele)

HaplotypeTree* init_haplotype_tree(int numSamples) {

    HaplotypeTree* tree = (HaplotypeTree*) calloc(1, sizeof(HaplotypeTree));

    tree -> numSamples = numSamples;
    tree -> leftHaplotype = (int*) calloc(numSamples, sizeof(int));
    tree -> rightHaplotype = (int*) calloc(numSamples, sizeof(int));

    tree -> labelMap = kh_init(32);

    tree -> numNodes = 1;
    tree -> prevNumGroups = 1;
    tree -> prevNumLeavesPerGroup = 1;

    return tree;

}

void add_locus(HaplotypeTree* tree, int numAlleles, GENOTYPE* genotypes, bool collapseMissingGenotypes) {

    for (int i = 0; i < tree -> numSamples; i++) {
        if (tree -> numNodes == 1) {
            tree -> leftHaplotype[i] = LEFT_ALLELE(genotypes[i]);
            tree -> rightHaplotype[i] = RIGHT_ALLELE(genotypes[i]);
        } else if (collapseMissingGenotypes && (tree -> leftHaplotype[i] == tree -> prevNumLeavesPerGroup || tree -> rightHaplotype[i] == tree -> prevNumLeavesPerGroup)) {
            tree -> leftHaplotype[i] = numAlleles + 1;
            tree -> rightHaplotype[i] = numAlleles + 1;
        } else {
            tree -> leftHaplotype[i] = CALCULATE_NEW_LABEL(tree -> leftHaplotype[i], tree -> numNodes, tree -> prevNumGroups, tree -> prevNumLeavesPerGroup, numAlleles, LEFT_ALLELE(genotypes[i]));
            tree -> rightHaplotype[i] = CALCULATE_NEW_LABEL(tree -> rightHaplotype[i], tree -> numNodes, tree -> prevNumGroups, tree -> prevNumLeavesPerGroup, numAlleles, RIGHT_ALLELE(genotypes[i])); 
        }
    }
    
    tree -> numNodes += (tree -> prevNumGroups) * (tree -> prevNumLeavesPerGroup) * (numAlleles + 1);
    tree -> prevNumGroups = tree -> prevNumLeavesPerGroup;
    tree -> prevNumLeavesPerGroup = (numAlleles + 1);

    if (tree -> numNodes >= MAX_NUM_NODES)
        relabel_haplotypes(tree);

}

void relabel_haplotypes(HaplotypeTree* tree) {

}

void destroy_haplotype_tree(HaplotypeTree* tree) {

    free(tree -> leftHaplotype);
    free(tree -> rightHaplotype);
    kh_destroy(32, tree -> labelMap);
    free(tree);

}

int main() {

}