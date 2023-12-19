
#ifndef _VCF_GENOTYPE_PARSER_
#define _VCF_GENOTYPE_PARSER_

#include <stdlib.h>

#include <stdbool.h>

#include "../klib/zlib.h"

#include "../klib/kstring.h"

#include "../klib/kseq.h"

#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

typedef char GENOTYPE;

typedef struct {
    kstring_t* file_name;
    gzFile file;
    kstream_t* stream;
    kstring_t* buffer;
    bool isEOF;

    int num_samples;
    kstring_t* sample_names;

    kstring_t* nextChromosome;
    int nextPosition;
    int nextNumAlleles;
    GENOTYPE* nextGenotypes;
} VCFGenotypeParser;

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name);

void get_next_locus(VCFGenotypeParser* parser, kstring_t* chromosome, int* position, int* numOfAlleles, GENOTYPE** genotypes);

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser);

static inline GENOTYPE parse_genotype(char* start, int numAlleles) {
    GENOTYPE genotype = (GENOTYPE) ((numAlleles + 1) << 4) | (numAlleles + 1);
    char* next = start + 1;
    if (start[0] != '.')
        genotype = (strtol(start, &next, 10) << 4) | (numAlleles + 1);
    if ((next[0] == '|' || next[0] == '/') && next[1] != '.')
        genotype = (genotype & 0xF0) | strtol(next + 1, (char**) NULL, 10);
    return genotype;
}

#endif