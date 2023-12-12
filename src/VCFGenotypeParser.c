
#include "VCFGenotypeParser.h"

#include <stdio.h>

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name) {

    if (!bgzf_is_bgzf(file_name))
        return NULL;

    BGZF* file = bgzf_open(file_name, "r");
    
    kstring_t* buffer = (kstring_t*) malloc(sizeof(kstring_t));
    do {
        bgzf_getline(file, '\n', buffer);
    } while (strncmp(ks_str(buffer), "#C", 2) != 0);
    
    int num_samples = 0;
    for (int i = 0; i < buffer -> l; i++)
        if (buffer -> s[i] == '\t')
            num_samples++;
    num_samples -= 8;

    kstring_t* sample_names = (kstring_t*) malloc(num_samples * sizeof(kstring_t));
    int num_tabs = 0, prev_index;
    for (int i = 0; i <= buffer -> l; i++) {
        if (i == buffer -> l || buffer -> s[i] == '\t') {
            if (num_tabs > 7) {
                kputsn(buffer -> s + i + 1, i - prev_index, &sample_names[num_tabs - 8]);
            }
            prev_index = i;
            num_tabs++;
        }
    }

    VCFGenotypeParser* parser = (VCFGenotypeParser*) malloc(sizeof(VCFGenotypeParser));
    parser -> file_name = (kstring_t*) malloc(sizeof(kstring_t));
    kputs(file_name, parser -> file_name);
    parser -> file = file;
    parser -> num_samples = num_samples;
    parser -> sample_names = sample_names;
    parser -> buffer = buffer;
    parser -> nextChromosome = (kstring_t*) malloc(sizeof(kstring_t));
    parser -> nextGenotypes = (GENOTYPE*) malloc(num_samples * sizeof(GENOTYPE));

    return parser;

}

bool get_next_locus(VCFGenotypeParser* parser, kstring_t* chromosome, int* position, int* numOfAlleles, GENOTYPE* genotypes) {
    return true;
}

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser) {
    free(parser -> file_name -> s);
    free(parser -> file_name);
    for (int i = 0; i < parser -> num_samples; i++) {
        free(parser -> sample_names[i].s);
    }
    free(parser -> sample_names);
    free(parser -> buffer -> s);
    free(parser -> buffer);
    free(parser -> nextChromosome -> s);
    free(parser -> nextChromosome);
    free(parser -> nextGenotypes);
}

int main() {
    init_vcf_genotype_parser("sample1.vcf.gz");
}