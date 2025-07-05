
// File: VCFLocusParser.c
// Date: 6 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parse a VCF file for samples' genotypes at each record.

#include "VCFLocusParser.h"

// I am not in love with how I wrote this, but it is sufficient for now.
bool seek(VCFLocusParser_t* parser) {

    int dret, numTabs, prevIndex, numAlleles;
    Locus l;
    double maf, afMissing, afMax, af;

    // We break out of the infinite loop until EOF is encountered or a 
    //  record that satisfies all the filters is encountered.
    while (true) {
        
        // Get next line.
        ks_getuntil(parser -> stream, '\n', parser -> buffer, &dret);

        // If EOF or nothing was read in (for safety), set EOF flag and return.
        if (isEOF(parser))
            return true;

        // This is alittle clunky, but I think it is faster than splitting on '\t'.
        numTabs = 0, prevIndex = 0, numAlleles = 2;
        for (int i = 0; i <= parser -> buffer -> l; i++) {
            // If end of the line is reached or a tab was encountered.
            if (i == parser -> buffer -> l || parser -> buffer -> s[i] == '\t') {
                if (numTabs == 0) {
                    // The first field in a record is the chromosome name.
                    if (parser -> nextChrom != NULL)
                        free(parser -> nextChrom);
                    parser -> nextChrom = strndup(parser -> buffer -> s, i);
                } else if (numTabs == 1) {
                    // The second field is the position on the chromosome.
                    parser -> nextCoord = (int) strtol(parser -> buffer -> s + prevIndex + 1, (char**) NULL, 10);
                } else if (numTabs == 4) {
                    // The fifth field holds the ALT alleles. Each record has at least two alleles.
                    //  Each additional allele is appended withparser -> buffer -> s a ','. For each ',' encountered,
                    //  increment the number of alleles in the record.
                    for (int j = prevIndex + 1; parser -> buffer -> s[j] != '\t'; j++)
                        if (parser -> buffer -> s[j] == ',')
                            numAlleles++;
                } else if (numTabs > 8) {
                    // The ninth field and on holds the genotypes of the samples.
                    l = parse_locus(parser -> buffer -> s + prevIndex + 1, numAlleles);
                    // Increment each allele's count.
                    parser -> alleleCounts[LEFT_ALLELE(l)]++;
                    parser -> alleleCounts[RIGHT_ALLELE(l)]++;
                    // Set the sample's corresponding genotype.
                    parser -> nextLocus[numTabs - 9] = l;
                }
                prevIndex = i;
                numTabs++;
            }
        }

        // Set the number of alleles at the locus.
        parser -> nextNumAlleles = numAlleles;
        
        // Calculate the number of missing genotypes present.
        afMissing = parser -> alleleCounts[numAlleles] / (2.0 * parser -> numSamples);
        parser -> alleleCounts[numAlleles] = 0;
        // Iterate through the allele counts to get the MAF and the most frequent allele.
        maf = 1;
        afMax = 0;
        for (int i = 0; i < numAlleles; i++) {
            af = parser -> alleleCounts[i] / (2.0 * parser -> numSamples - parser -> alleleCounts[numAlleles]);
            if (af < maf)
                maf = af;
            if (af > afMax)
                afMax = af;
            parser -> alleleCounts[i] = 0;
        }

        // Test that all thresholds are met.
        if (numAlleles > MAX_NUM_ALLELES)
            continue;
        else if (parser -> dropMonomorphicSites && afMax == 1)
            continue;
        else if (numAlleles == 2 && maf < parser -> maf)
            continue;
        else if (afMissing > parser -> afMissing)
            continue;
        else 
            break;

    }

    return false;

}

VCFLocusParser_t* init_vcf_locus_parser(char* fileName, double maf, double afMissing, bool dropMonomorphicSites) {

    // Open the GZ file.
    gzFile file = gzopen(fileName, "r");

    // If file does not exist or is not compressed using gzip, return NULL.
    int errnum;
    gzerror(file, &errnum);
    if (errnum != Z_OK) {
        return NULL;
    }
    
    // Initialize the file stream.
    kstream_t* stream = ks_init(file);
    
    // Initialize the buffer to read in from the stream.
    kstring_t* buffer = calloc(1, sizeof(kstring_t));
    
    // Parse all the meta data in the VCF file.
    int dret;
    do {
        ks_getuntil(stream, '\n', buffer, &dret);
    } while (strncmp(buffer -> s, "#C", 2) != 0);
    
    // Count the number of samples in the header of the VCF file.
    int numSamples = 0;
    for (int i = 0; i < buffer -> l; i++)
        if (buffer -> s[i] == '\t')
            numSamples++;
    numSamples -= 8;
    
    // Allocate the array to hold the sample names.
    char** sampleNames = (char**) calloc(numSamples, sizeof(char*));
    // Read in the sample names.
    int numTabs = 0, prevIndex;
    for (int i = 0; i <= buffer -> l; i++) {
        if (i == buffer -> l || buffer -> s[i] == '\t') {
            if (numTabs > 8)
                sampleNames[numTabs - 9] = strndup(buffer -> s + prevIndex + 1, i - prevIndex - 1);
            prevIndex = i;
            numTabs++;
        }
    }
    
    // Allocate our structure and the memory for its fields.
    VCFLocusParser_t* parser = (VCFLocusParser_t*) calloc(1, sizeof(VCFLocusParser_t));
    parser -> fileName = strdup(fileName);
    parser -> file = file;
    parser -> stream = stream;
    parser -> numSamples = numSamples;
    parser -> sampleNames = sampleNames;
    parser -> buffer = buffer;
    parser -> nextChrom = NULL;
    parser -> nextLocus = (Locus*) calloc(numSamples, sizeof(Locus));

    // Set fields from arguments.
    parser -> maf = maf;
    parser -> afMissing = afMissing;
    parser -> dropMonomorphicSites = dropMonomorphicSites;
    parser -> alleleCounts = calloc(MAX_NUM_ALLELES, sizeof(int));
    for (int i = 0; i < MAX_NUM_ALLELES; i++)
        parser -> alleleCounts[i] = 0;

    // Prime the first record.
    seek(parser);
    
    // Return created parser.
    return parser;
}

bool get_next_locus(VCFLocusParser_t* parser, char** chrom, int* coord, int* numOfAlleles, Locus** locus) {
    
    // Move primed read into arguments.
    if (*chrom != NULL)
        free(*chrom);
    *chrom = strdup(parser -> nextChrom);
    *coord = parser -> nextCoord;
    *numOfAlleles = parser -> nextNumAlleles;
    // We just swap array pointers, which is better than copying each element individually.
    Locus* temp = *locus;
    *locus = parser -> nextLocus;
    parser -> nextLocus = temp;
    
    // Prime the next read.
    return seek(parser);
}

bool isEOF(VCFLocusParser_t* parser) {
    return ks_eof(parser -> stream) || parser -> buffer -> l == 0;
}

void destroy_vcf_locus_parser(VCFLocusParser_t* parser) {
    if (parser == NULL)
        return;
    // Destroy the buffer.
    if (parser -> buffer != NULL) {
        if (parser -> buffer -> s != NULL)
            free(parser -> buffer -> s);
        free(parser -> buffer);
    }
    // Close the file being read in.
    gzclose(parser -> file);
    // Destroy the stream.
    ks_destroy(parser -> stream);
    // Free the file name.
    free(parser -> fileName);
    // Destroy the sample names.
    for (int i = 0; i < parser -> numSamples; i++)
        if (parser -> sampleNames[i] != NULL)
            free(parser -> sampleNames[i]);
    free(parser -> sampleNames);
    // Destroy the next chromosome string.
    if (parser -> nextChrom != NULL)
        free(parser -> nextChrom);
    // Destroy the locus array.
    free(parser -> nextLocus);
    // Free allele counts.
    free(parser -> alleleCounts);
    // Destroy the parser.
    free(parser);
}