
// File: VCFGenotypeParser.h
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Parse the genotypes for each record in a VCF file.

#ifndef _VCF_GENOTYPE_PARSER_
#define _VCF_GENOTYPE_PARSER_

#include <stdlib.h>

#include <stdbool.h>

// Supports GZIP VCF files.
#include "../klib/zlib.h"

#include "../klib/kstring.h"

// We use the klib wrapper to read in streams.
#include "../klib/kseq.h"

// Initiate the klib stream.
#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

// A sample's genotype is encoded in a byte.
//  Therefore, there is a maximum of 15 possible
//  alleles and the missing allele at each locus.
typedef char GENOTYPE;

// Our parser structure.
typedef struct {
    // The name of the VCF file.
    kstring_t* file_name;
    // Our GZ file structure.
    gzFile file;
    // The stream we will read from the gzFile.
    kstream_t* stream;
    // The dynamic buffer used by kseq.
    kstring_t* buffer;
    // Flag set when EOF.
    bool isEOF;

    // The number of samples in the VCF file.
    int num_samples;
    // The names of the samples.
    kstring_t* sample_names;

    // By adding a "peak" mechanism to the VCF file, such as we can do in streams,
    //  many algorithms are simplified. The stream is pointing the record after
    //  these entries.
    kstring_t* nextChromosome;
    int nextPosition;
    int nextNumAlleles;
    GENOTYPE* nextGenotypes;
} VCFGenotypeParser;

// Creates a VCFGenotypeParser.
// Accepts:
//  char* file_name -> The name of the file to read in.
// Returns:
//  The created parser or NULL if file does not exist.
VCFGenotypeParser* init_vcf_genotype_parser(char* file_name);

// Get the next record from a parser.
// Accepts:
//  VCFGenotypeParser* parser -> A pointer to a parser.
//  kstring_t* chromosome -> Sets the chromosome of the record.
//  int* position -> Sets the position of the record.
//  int* numOfAlleles -> Sets the number of alleles at that locus.
//  GENOTYPE** genotypes -> Fills an array of samples' genotypes. Is set (swapped)
//                              with the nextGenotypes array in VCFGenotypeParser.
// Returns: 
//  void. Pointers are left unchanged when isEOF.
void get_next_locus(VCFGenotypeParser* parser, kstring_t* chromosome, int* position, int* numOfAlleles, GENOTYPE** genotypes);

// Deallocate all the memory occupied by the VCFGenotypeParser.
// Accepts:
//  VCFGenotypeParser* parser -> The parser to destroy.
// Returns:
//  void.
void destroy_vcf_genotype_parser(VCFGenotypeParser* parser);

// Encodes the genotype of a sample into a byte.
//  The function that is called the most. Is it fast enough?
// Accepts:
//  char* start -> A pointer to the beginning of a genotype in a VCF record.
//  int numAlleles -> The number of possible alleles at the record.
// Returns:
//  GENOTYPE, The encoded genotype of the sample. 
static inline GENOTYPE parse_genotype(char* start, int numAlleles) {
    // There are numAlleles at the locus labeled 0 ... numAlleles - 1.
    //  The numAllele denotes the missing allele. Set initial genotype to
    //  two missing alleles.
    GENOTYPE genotype = (GENOTYPE) (numAlleles << 4) | numAlleles;
    char* next = start + 1;
    // If the left allele is not missing, then parse integer and set left genotype.
    if (start[0] != '.')
        genotype = (strtol(start, &next, 10) << 4) | numAlleles;
    // If there is a second, non-missing genotype, parse and set right genotype.
    if ((next[0] == '|' || next[0] == '/') && next[1] != '.')
        genotype = (genotype & 0xF0) | strtol(next + 1, (char**) NULL, 10);
    return genotype;
}

#endif