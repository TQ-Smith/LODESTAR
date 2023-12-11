
#ifndef _VCF_GENOTYPE_PARSER_
#define _VCF_GENOTYPE_PARSER_

#include <stdlib.h>

#include <stdbool.h>

#include "../klib/kstring.h"

#include "../klib/bgzf.h"

#include "../klib/kdq.h"

KDQ_INIT(kstring_t)

typedef char GENOTYPE;

typedef struct {
    kstring_t file_name;
    BGZF* file;
    int num_samples;
    kdq_t(kstring_t) *sample_names;

    GENOTYPE* nextGenotypes;
    int nextPosition;
    kstring_t nextChromosome;
} VCFGenotypeParser;

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name);

bool getNextLocus(VCFGenotypeParser* parser, char* chromosome, int* position, GENOTYPE* genotypes);

char** getSampleNames(VCFGenotypeParser* parser);

destroy_vcf_genotype_parser(VCFGenotypeParser* parser);

#endif