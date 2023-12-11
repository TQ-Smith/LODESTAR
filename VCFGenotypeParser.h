
#ifndef _VCF_GENOTYPE_PARSER_
#define _VCF_GENOTYPE_PARSER_

#include <stdlib.h>

#include <string.h>

#include "bgzf.h"

#include "kdq.h"

KDQ_INIT(char*)

typedef char* GENOTYPE;

typedef struct VCFGenotypeParser {
    char* file_name;
    BGZF* file;
    int num_samples;
    kdq_t(char*) *sample_names;

    GENOTYPE* nextGenotypes;
    int nextPosition;
    char nextChromosome[15];
} VCFGenotypeParser;

#endif