
// File: VCFLocusParser.h
// Date: 6 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parse a VCF file for samples' genotypes at each record.

#ifndef _VCF_LOCUS_PARSER_H_
#define _VCF_LOCUS_PARSER_H_

#include "../lib/zlib.h"
#include "../lib/kstring.h"
#include "../lib/kseq.h"
#include "RegionSet.h"

// We use kseq as a stream to read in GZ files.
#define BUFFER_SIZE 4096
KSTREAM_INIT(gzFile, gzread, BUFFER_SIZE)

// We define a sample's genotype as a char.
//  The most significant 4-bits are the left allele, 
//  and the least significant 4-bits are the right allele.
//  NOTE: If a locus has N alleles labeled 0,..,N - 1, then the
//  missing allele is encoded as N. Thus, there is a maximum of 2^16
//  alleles at each locus. This can be easily modified. Since this number
//  is large, we artifically set the maximum number of alleles at a locus to 64.
typedef unsigned int Locus;
#define LEFT_ALLELE(a) (a >> 16)
#define RIGHT_ALLELE(a) (a & 0x0000FFFF)
#define MAX_NUM_ALLELES 64

// Our structure that represents a VCFLocusParser.
typedef struct {
    // The name of the VCF file.
    kstring_t* fileName;
    // The GZ file. Note: Seamlessly works with uncompressed files too.
    gzFile file;
    // The stream to read in the GZ file.
    kstream_t* stream;
    // A buffer to hold a line of the VCF parser from the stream.
    kstring_t* buffer;
    bool isEOF;

    // Number of samples in the VCF file.
    int numSamples;
    // An array of kstring_t to hold the names of each sample.
    kstring_t** sampleNames;

    // Fields to help us omitt unwanted record.
    //  If the user wishes, they can specify a set of regions to include/exclude.
    RegionSet_t* set;
    // Flag to drop monomorphic sites.
    bool dropMonomorphicSites;
    // Minor-allele-frequency threshold.
    double maf;
    // Missing genotype threshold.
    double afMissing;
    // An array to hold the number of each allele at a locus.
    //  Used in maf calculation.
    int alleleCounts[MAX_NUM_ALLELES];

    // For convenience, we implement a priming read/peak operation.
    //  This allows us to easily test if the next record belongs to a different chromosome.
    kstring_t* nextChrom;
    unsigned int nextCoord;
    int nextNumAlleles;
    // Array that holds the genotypes for each of the samples.
    Locus* nextLocus;
} VCFLocusParser_t;

// Creates a VCFLocusParser_t structure.
// Accepts:
//  char* fileName -> The name of the VCF file.
//  kstring_t* regions -> Regions to include/exclude when reading in the VCF file.
//                          See RegionSet.h on how to specify string.
//                          All regions are parsed if NULL.
//  bool takeComplement -> When false, regions in regions parameter are included during parsing.
//                          When true, regions in the regions parameter are excluded during parsing.
//                          Not used if regions is NULL.
//  double maf -> Records with a minor-allele-frequency < maf are dropped.
//  double afMissing -> Record with a proportion of missing genotypes >= afMissing are dropped.
//  bool dropMonomorphicSites -> If set, records with monomorphic genotypes are dropped.
// Returns: VCFLocusParser_t*, If fileName does not correspond to a file or regions could not be parsed, return NULL.
//              Otherwise, return the created structure.
VCFLocusParser_t* init_vcf_locus_parser(char* fileName, kstring_t* regions, bool takeComplement, double maf, double afMissing, bool dropMonomorphicSites);

// Get the next locus from a parser.
// Accepts:
//  VCFLocusParser_t* parser -> The parser structure to read the VCF file.
//  kstring_t* chrom -> Sets the chromosome of the read-in record.
//  unsigned int* coord -> Sets the position of the read-in record.
//  int* numOfAlleles -> Sets the number of alleles at the read-in record.
//  Locus** genos -> Sets the array of the samples' genotypes at the read-in record.
// Returns: void. Note: when EOF is set, passed arguments are unchanged.
void get_next_locus(VCFLocusParser_t* parser, kstring_t* chrom, unsigned int* coord, int* numOfAlleles, Locus** genos);

// Free all the memory allocated to a VCFLocusParser_t.
// Accepts:
//  VCFLocusParser_t* parser -> The parser to free.
// Returns: void.
void destroy_vcf_locus_parser(VCFLocusParser_t* parser);

// Encode a string genotype from the VCF file into a Locus.
//  This function will be called a lot.
// Accepts:
//  char* start -> The genotype to parse.
//  int numAlleles -> The number alleles at the given record.
// Returns: Locus, The encoded genotype.
static inline Locus parse_locus(char* start, int numAlleles) {
    // By default, both genotypes at a locus are missing.
    Locus locus = (numAlleles << 16) | numAlleles;
    char* next = start + 1;
    // If the left allele is not missing, then parse integer and set left genotype.
    if (start[0] != '.')
        locus = (strtol(start, &next, 10) << 16) | numAlleles;
    // If there is a second, non-missing genotype, parse and set right genotype.
    if ((next[0] == '|' || next[0] == '/') && next[1] != '.')
        locus = (locus & 0xFFFF0000) | strtol(next + 1, (char**) NULL, 10);
    return locus;
}

#endif