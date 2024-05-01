
// File: 
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Parse the genotypes for each record in a VCF file.

#ifndef _VCF_LOCUS_PARSER_H_
#define _VCF_LOCUS_PARSER_H_

#include <stdlib.h>

#include <stdbool.h>

#include "../lib/zlib.h"

#include "../lib/kstring.h"

#include "../lib/kseq.h"

#include "RegionFilter.h"

#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

typedef char Locus;

#define LEFT_ALLELE(a) (a >> 4)
#define RIGHT_ALLELE(a) (a & 0x0F)

typedef struct {
    kstring_t* fileName;
    gzFile file;
    kstream_t* stream;
    kstring_t* buffer;
    bool isEOF;

    int numSamples;
    kstring_t* sampleNames;

    RegionFilter* filter;
    double maf;
    double afMissing;
    int alleleCounts[16];

    kstring_t* nextChrom;
    unsigned int nextPos;
    int nextNumAlleles;
    Locus* nextLocus;
} VCFLocusParser;

VCFLocusParser* init_vcf_locus_parser(char* fileName, RegionFilter* filter, double maf, double afMissing);

void get_next_locus(VCFLocusParser* parser, kstring_t* chrom, unsigned int* pos, int* numOfAlleles, Locus** genos);

void destroy_vcf_locus_parser(VCFLocusParser* parser);

static inline Locus parse_locus(char* start, int numAlleles) {
    Locus locus = (Locus) (numAlleles << 4) | numAlleles;
    char* next = start + 1;
    // If the left allele is not missing, then parse integer and set left genotype.
    if (start[0] != '.')
        locus = (strtol(start, &next, 10) << 4) | numAlleles;
    // If there is a second, non-missing genotype, parse and set right genotype.
    if ((next[0] == '|' || next[0] == '/') && next[1] != '.')
        locus = (locus & 0xF0) | strtol(next + 1, (char**) NULL, 10);
    return locus;
}

#endif