
#ifndef _VCF_GENOTYPE_PARSER_
#define _VCF_GENOTYPE_PARSER_

#include <stdlib.h>

#include <stdbool.h>

#include "../klib/kstring.h"

#include "../klib/bgzf.h"

#include "../klib/kstring.h"

typedef char GENOTYPE;

typedef struct {
    kstring_t* file_name;
    BGZF* file;
    int num_samples;
    kstring_t* sample_names;

    kstring_t* buffer;

    kstring_t* nextChromosome;
    int nextPosition;
    int nextNumAlleles;
    GENOTYPE* nextGenotypes;
} VCFGenotypeParser;

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name);

bool get_next_locus(VCFGenotypeParser* parser, kstring_t* chromosome, int* position, int* numOfAlleles, GENOTYPE* genotypes);

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser);

#endif