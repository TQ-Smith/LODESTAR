
// File: VCFLocusParser.c
// Date: 6 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Parse a VCF file for samples' genotypes at each record.

#include "VCFLocusParser.h"

VCFLocusParser_t* init_vcf_locus_parser(char* fileName, kstring_t* regions, bool takeComplement, double maf, double afMissing, bool dropMonomorphicSites) {

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
    kstring_t* buffer = init_kstring(NULL);
    
    // Parse all the meta data in the VCF file.
    int dret;
    do {
        ks_getuntil(stream, '\n', buffer, &dret);
    } while (strncmp(ks_str(buffer), "#C", 2) != 0);
    
    // Count the number of samples in the header of the VCF file.
    int numSamples = 0;
    for (int i = 0; i < ks_len(buffer); i++)
        if (buffer -> s[i] == '\t')
            numSamples++;
    numSamples -= 8;
    
    // Allocate the array to hold the sample names.
    kstring_t** sampleNames = (kstring_t**) calloc(numSamples, sizeof(kstring_t*));
    // Read in the sample names.
    int numTabs = 0, prevIndex;
    for (int i = 0; i <= ks_len(buffer); i++) {
        if (i == ks_len(buffer) || ks_str(buffer)[i] == '\t') {
            if (numTabs > 8) {
                sampleNames[numTabs - 9] = init_kstring(NULL);
                kputsn(ks_str(buffer) + prevIndex + 1, i - prevIndex - 1, sampleNames[numTabs - 9]);
            }
            prevIndex = i;
            numTabs++;
        }
    }
    
    // Allocate our structure and the memory for its fields.
    VCFLocusParser_t* parser = (VCFLocusParser_t*) calloc(1, sizeof(VCFLocusParser_t));
    parser -> fileName = init_kstring(fileName);
    parser -> file = file;
    parser -> stream = stream;
    parser -> numSamples = numSamples;
    parser -> sampleNames = sampleNames;
    parser -> buffer = buffer;
    parser -> isEOF = false;
    parser -> nextChrom = init_kstring(NULL);
    parser -> nextLocus = (Locus*) calloc(numSamples, sizeof(Locus));
    if (regions == NULL)
        parser -> set = NULL;
    else {
        // Try parsing the regions argument.
        parser -> set = init_region_set(regions, takeComplement);
        // If error, then free memory already allocated to parser, and return NULL for error.
        if (parser -> set == NULL) {
            destroy_vcf_locus_parser(parser);
            return NULL;
        }
    }

    // Set fields from arguments.
    parser -> maf = maf;
    parser -> afMissing = afMissing;
    parser -> dropMonomorphicSites = dropMonomorphicSites;
    for (int i = 0; i < 16; i++)
        parser -> alleleCounts[i] = 0;

    // Prime the first record.
    get_next_locus(parser, parser -> nextChrom, &(parser -> nextCoord), &(parser -> nextNumAlleles), &(parser -> nextLocus));
    
    // Return created parser.
    return parser;
}

// I am not in love with how I wrote this, but it is sufficient for now.
void seek(VCFLocusParser_t* parser) {

    int dret, numTabs, prevIndex, numAlleles;
    bool isInSet;
    Locus l;
    double maf, afMissing, afMax, af;

    // We break out of the infinite loop until EOF is encountered or a 
    //  record that satisfies all the filters is encountered.
    while (true) {

        // Get next line.
        ks_getuntil(parser -> stream, '\n', parser -> buffer, &dret);

        // If EOF or nothing was read in (for safety), set EOF flag and return.
        if (ks_eof(parser -> stream) || ks_len(parser -> buffer) == 0) {
            parser -> isEOF = true;
            return;
        }

        // This is alittle clunky, but I think it is faster than splitting on '\t'.
        numTabs = 0, prevIndex = 0, numAlleles = 2;
        isInSet = true;
        for (int i = 0; i <= ks_len(parser -> buffer) && isInSet; i++) {
            // If end of the line is reached or a tab was encountered.
            if (i == ks_len(parser -> buffer) || ks_str(parser -> buffer)[i] == '\t') {
                if (numTabs == 0) {
                    // The first field in a record is the chromosome name.
                    ks_overwriten(ks_str(parser -> buffer), i, parser -> nextChrom);
                } else if (numTabs == 1) {
                    // The second field is the position on the chromosome.
                    parser -> nextCoord = (int) strtol(ks_str(parser -> buffer) + prevIndex + 1, (char**) NULL, 10);
                    // If a RegionSet was specified and the locus is not in the set, set flag and stop parsing record.
                    if (parser -> set != NULL && !query_locus(parser -> set, parser -> nextChrom, (unsigned int) parser -> nextCoord))
                        isInSet = false;
                } else if (numTabs == 4) {
                    // The fifth field holds the ALT alleles. Each record has at least two alleles.
                    //  Each additional allele is appended with a ','. For each ',' encountered,
                    //  increment the number of alleles in the record.
                    for (int j = prevIndex + 1; ks_str(parser -> buffer)[j] != '\t'; j++)
                        if (ks_str(parser -> buffer)[j] == ',')
                            numAlleles++;
                    // If the number of alleles exceeds the maximum, drop record.
                    if (numAlleles > MAX_NUM_ALLELES)
                        isInSet = false;
                } else if (numTabs > 8) {
                    // The ninth field and on holds the genotypes of the samples.
                    l = parse_locus(ks_str(parser -> buffer) + prevIndex + 1, numAlleles);
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

        // If the record is not in the set, skip record.
        if (!isInSet)
            continue;
        
        // Calculate the number of missing genotypes present.
        afMissing = parser -> alleleCounts[numAlleles] / (2.0 * parser -> numSamples);
        parser -> alleleCounts[numAlleles] = 0;
        // Iterate through the allele counts to get the MAF and the most frequent allele.
        maf = 1;
        afMax = 0;
        for (int i = 0; i < numAlleles; i++) {
            af = parser -> alleleCounts[i] / (2.0 * parser -> numSamples);
            if (af < maf)
                maf = af;
            if (af > afMax)
                afMax = af;
            parser -> alleleCounts[i] = 0;
        }

        // Test that all thresholds are met.
        if (parser -> dropMonomorphicSites && afMax == 1)
            continue;
        else if (numAlleles == 2 && maf < parser -> maf)
            continue;
        else if (afMissing > parser -> afMissing)
            continue;
        else 
            return;

    }

}

void get_next_locus(VCFLocusParser_t* parser, kstring_t* chrom, unsigned int* coord, int* numOfAlleles, Locus** locus) {
    // If parser does not exist or EOF, leave arguments unchanged and return.
    if (parser == NULL || parser -> isEOF)
        return;
    
    // Move the primed read into the arguments.
    if (parser -> nextChrom != chrom)
        ks_overwrite(ks_str(parser -> nextChrom), chrom);
    *coord = parser -> nextCoord;
    *numOfAlleles = parser -> nextNumAlleles;
    // We just swap array pointers, which is better than copying each element individually.
    Locus* temp = *locus;
    *locus = parser -> nextLocus;
    parser -> nextLocus = temp;

    // Prime the next read.
    seek(parser);
}


void destroy_vcf_locus_parser(VCFLocusParser_t* parser) {
    if (parser == NULL)
        return;
    // Close the file being read in.
    gzclose(parser -> file);
    // Destroy the stream.
    ks_destroy(parser -> stream);
    // Free the file name.
    destroy_kstring(parser -> fileName);
    // If region set was allocated, destroy set.
    if (parser -> set != NULL)
        destroy_region_set(parser -> set);
    // Destroy the sample names.
    for (int i = 0; i < parser -> numSamples; i++)
        destroy_kstring(parser -> sampleNames[i]);
    free(parser -> sampleNames);
    // Destroy the buffer.
    destroy_kstring(parser -> buffer);
    // Destroy the next chromosome string.
    destroy_kstring(parser -> nextChrom);
    // Destroy the locus array.
    free(parser -> nextLocus);
    // Destroy the parser.
    free(parser);
}

// Used to test the parser.
/*
int main() {
    
    kstring_t* intervals = (kstring_t*) calloc(1, sizeof(kstring_t));
    ks_overwrite("3", intervals);
    VCFLocusParser_t* parser = init_vcf_locus_parser("./data/vcf_parser_test.vcf.gz", NULL, false, 0.4, 0.1, true);

    printf("There are %d samples with the following names:\n", parser -> numSamples);
    for (int i = 0; i < parser -> numSamples; i++)
        printf("%s\n", ks_str(&(parser -> sampleNames[i])));
    
    kstring_t* chromosome = (kstring_t*) calloc(1, sizeof(kstring_t));
    unsigned int position;
    int numOfAlleles;
    Locus* locus = (Locus*) calloc(parser -> numSamples, sizeof(Locus));

    while(!parser -> isEOF) {

        get_next_locus(parser, chromosome, &position, &numOfAlleles, &locus);
        
        printf("\n");
        printf("Chromosome: %s\n", ks_str(chromosome));
        printf("Position: %d\n", position);
        printf("Num Alleles: %d\n", numOfAlleles);
        printf("Genotypes:\n");
        for (int i = 0; i < parser -> numSamples; i++)
            printf("%x\n", locus[i]);
        printf("\n");
        
    }

    destroy_vcf_locus_parser(parser);
    free(locus);
    free(ks_str(chromosome)); free(chromosome);
    free(ks_str(intervals)); free(intervals);
    
}
*/