
#include "VCFGenotypeParser.h"

#include <stdio.h>
#include <assert.h>

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name) {

    if (!bgzf_is_bgzf(file_name))
        return NULL;

    BGZF* file = bgzf_open(file_name, "r");
    
    kstring_t* buffer = (kstring_t*) malloc(sizeof(kstring_t));
    do {
        bgzf_getline(file, '\n', buffer);
    } while (strncmp(ks_str(buffer), "#C", 2) != 0);
    
    int num_samples = 0;
    for (int i = 0; i < ks_len(buffer); i++)
        if (ks_str(buffer)[i] == '\t')
            num_samples++;
    num_samples -= 8;

    char** sample_names = (char**) malloc(num_samples * sizeof(char*));
    int prev_index = -1;
    for (int i = 0; i < ks_len(buffer); i++) {
    }

    VCFGenotypeParser* parser = (VCFGenotypeParser*) malloc(sizeof(VCFGenotypeParser));
    parser -> file_name = file_name;
    parser -> file = file;
    parser -> num_samples = num_samples;
    parser -> sample_names = sample_names;

    return parser;

}

bool getNextLocus(VCFGenotypeParser* parser, char* chromosome, int* position, GENOTYPE* genotypes) {
    return true;
}

char** getSampleNames(VCFGenotypeParser* parser) {
    return NULL;
}

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser) {

}

int main() {
    init_vcf_genotype_parser("sample1.vcf.gz");
}