
#ifndef _VCF_GENOTYPE_PARSER_
#define _VCF_GENOTYPE_PARSER_

#include <stdlib.h>

#include <stdbool.h>

#include "../klib/kstring.h"

#include "../klib/bgzf.h"

#include "../klib/kstring.h"

typedef char GENOTYPE;

typedef struct {
    char* file_name;
    BGZF* file;
    int num_samples;
    char** sample_names;

    GENOTYPE* nextGenotypes;
    int nextPosition;
    char* nextChromosome;
} VCFGenotypeParser;

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name);

bool getNextLocus(VCFGenotypeParser* parser, char* chromosome, int* position, GENOTYPE* genotypes);

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser);

#endif