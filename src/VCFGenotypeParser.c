
#include "VCFGenotypeParser.h"

#include <stdio.h>

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name) {
    
    if (!bgzf_is_bgzf(file_name)) {
        return NULL;
    }

    BGZF* file = bgzf_open(file_name, "r");
    
    kstring_t* buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
    do {
        bgzf_getline(file, '\n', buffer);
    } while (strncmp(ks_str(buffer), "#C", 2) != 0);
    
    int num_samples = 0;
    for (int i = 0; i < ks_len(buffer); i++)
        if (buffer -> s[i] == '\t')
            num_samples++;
    num_samples -= 8;

    kstring_t* sample_names = (kstring_t*) calloc(num_samples, sizeof(kstring_t));
    int num_tabs = 0, prev_index;
    for (int i = 0; i <= ks_len(buffer); i++)
        if (ks_str(buffer)[i] == '\t') {
            if (num_tabs > 7)
                kputsn(ks_str(buffer) + i + 1, i - prev_index, &sample_names[num_tabs - 8]);
            prev_index = i;
            num_tabs++;
        }

    VCFGenotypeParser* parser = (VCFGenotypeParser*) malloc(sizeof(VCFGenotypeParser));
    parser -> file_name = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(file_name, parser -> file_name);
    parser -> file = file;
    parser -> num_samples = num_samples;
    parser -> sample_names = sample_names;
    parser -> buffer = buffer;
    parser -> isEOF = false;
    parser -> nextChromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    parser -> nextGenotypes = (GENOTYPE*) calloc(num_samples, sizeof(GENOTYPE));

    get_next_locus(parser, parser -> nextChromosome, &(parser -> nextPosition), &(parser -> nextNumAlleles), &(parser -> nextGenotypes));

    return parser;

}

void get_next_locus(VCFGenotypeParser* parser, kstring_t* chromosome, int* position, int* numOfAlleles, GENOTYPE** genotypes) {

    if (parser == NULL || parser -> isEOF)
        return;

    if (bgzf_check_EOF(parser -> file)) {
        parser -> isEOF = true;
        return;
    }

    kputs(ks_str(parser -> nextChromosome), chromosome);
    *position = parser -> nextPosition;
    *numOfAlleles = parser -> nextNumAlleles;
    GENOTYPE* temp = *genotypes;
    *genotypes = parser -> nextGenotypes;
    parser -> nextGenotypes = temp;

    bgzf_getline(parser -> file, '\n', parser -> buffer);

    int num_tabs = 0, prev_index, numAlleles = 2;
    for (int i = 0; i <= ks_len(parser -> buffer); i++) {
        if (ks_str(parser -> buffer)[i] == '\t') {
            if (num_tabs == 0)
                kputsn(ks_str(parser -> buffer) + i + 1, i - prev_index, parser -> nextChromosome);
            else if (num_tabs == 1)
                parser -> nextPosition = (int) strtol(ks_str(parser -> buffer) + i + 1, (char**) NULL, 10);
            else if (num_tabs == 4) {
                for (int j = i + 1; ks_str(parser -> buffer)[j] != '\t'; j++)
                    if (ks_str(parser -> buffer)[j] == ',')
                        numAlleles++;
            } else if (num_tabs > 7)
                parser -> nextGenotypes[num_tabs - 8] = parse_genotype(ks_str(parser -> buffer) + i + 1);
            prev_index = i;
            num_tabs++;
        }
    }

    *numOfAlleles = numAlleles;

}

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser) {
    if (parser == NULL)
        return;

    bgzf_close(parser -> file);
    free(ks_str(parser -> file_name)); free(parser -> file_name);
    for (int i = 0; i < parser -> num_samples; i++)
        free(parser -> sample_names[i].s);
    free(parser -> sample_names);
    free(ks_str(parser -> buffer)); free(parser -> buffer);
    free(ks_str(parser -> nextChromosome)); free(parser -> nextChromosome);
    free(parser -> nextGenotypes);
}

int main() {
    
    VCFGenotypeParser* parser = init_vcf_genotype_parser("vcf_parser_test.vcf.gz");

    printf("There are %d samples with the following names:\n", parser -> num_samples);
    for (int i = 0; i < parser -> num_samples; i++)
        printf("%s\n", ks_str(&(parser -> sample_names[i])));
    printf("%d", bgzf_check_EOF(parser -> file));
    /*
    kstring_t* chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    int position;
    int numOfAlleles;
    GENOTYPE* genotypes = (GENOTYPE*) calloc(parser -> num_samples, sizeof(GENOTYPE));

    while(true) {

        get_next_locus(parser, chromosome, &position, &numOfAlleles, &genotypes);

        if (parser -> isEOF)
            break;
        
        printf("\n");
        printf("Chromosome: %s\n", ks_str(chromosome));
        printf("Position: %d\n", position);
        printf("Num Alleles: %d\n", numOfAlleles);
        printf("Genotypes:\n");
        for (int i = 0; i < parser -> num_samples; i++)
            printf("%x\n", genotypes[i]);
        printf("\n");
        
    }

    destroy_vcf_genotype_parser(parser);
    free(genotypes);
    free(ks_str(chromosome)); free(chromosome);
    */
}