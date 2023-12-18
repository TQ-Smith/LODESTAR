
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
    bool isEOF;

    kstring_t* nextChromosome;
    int nextPosition;
    int nextNumAlleles;
    GENOTYPE* nextGenotypes;
} VCFGenotypeParser;

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name);

void get_next_locus(VCFGenotypeParser* parser, kstring_t* chromosome, int* position, int* numOfAlleles, GENOTYPE** genotypes);

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser);

static inline GENOTYPE parse_genotype(char* start) {
    GENOTYPE genotype = (GENOTYPE) 0xFF;
    char* next = NULL;
    if (start[0] != '.')
        genotype = ((char) strtol(start, &next, 10) << 4) | (char) 0x0F;
    if (next[1] != '.' && next[1] != '\t' && next[1] != ':')
        genotype &= (0xF0 | (char) strtol(next + 1, (char**) NULL, 10));
    return genotype;
}

#endif