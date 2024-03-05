
// File: VCFGenotypeParser.c
// Date: 18 Janurary 2024
// Author: TQ Smith
// Purpose: Parse the genotypes for each record in a VCF file.

#include <stdio.h>

#include "VCFGenotypeParser.h"

VCFGenotypeParser* init_vcf_genotype_parser(char* fileName) {

    // Try to open file.
    gzFile file = gzopen(fileName, "r");

    // Make sure the file is in GZ format.
    int errnum;
    gzerror(file, &errnum);
    if (errnum != Z_OK)
        return NULL;
    
    // Create stream.
    kstream_t* stream = ks_init(file);
    
    // Create buffer to read in VCF file record-by-record.
    kstring_t* buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
    int dret;

    // Swallow header lines.
    do {
        ks_getuntil(stream, '\n', buffer, &dret);
    } while (strncmp(ks_str(buffer), "#C", 2) != 0);
    
    // Count the number of samples.
    int numSamples = 0;
    for (int i = 0; i < ks_len(buffer); i++)
        if (buffer -> s[i] == '\t')
            numSamples++;
    numSamples -= 8;

    // Allocate array to hold sample names and fill array.
    kstring_t* sampleNames = (kstring_t*) calloc(numSamples, sizeof(kstring_t));
    int numTabs = 0, prevIndex;
    for (int i = 0; i <= ks_len(buffer); i++)
        if (ks_str(buffer)[i] == '\t') {
            if (numTabs > 7)
                kputsn(ks_str(buffer) + i + 1, i - prevIndex, &sampleNames[numTabs - 8]);
            prevIndex = i;
            numTabs++;
        }

    // Allocate all necessary parser memory and set values.
    VCFGenotypeParser* parser = (VCFGenotypeParser*) malloc(sizeof(VCFGenotypeParser));
    parser -> fileName = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(fileName, parser -> fileName);
    parser -> file = file;
    parser -> stream = stream;
    parser -> numSamples = numSamples;
    parser -> sampleNames = sampleNames;
    parser -> buffer = buffer;
    parser -> isEOF = false;
    parser -> nextChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    parser -> nextGenos = (GENOTYPE*) calloc(numSamples, sizeof(GENOTYPE));

    // Read in first locus to prime the read.
    get_next_locus(parser, parser -> nextChrom, &(parser -> nextPos), &(parser -> nextNumAlleles), &(parser -> nextGenos));

    // Return structure.
    return parser;

}

void get_next_locus(VCFGenotypeParser* parser, kstring_t* chrom, int* pos, int* numOfAlleles, GENOTYPE** genos) {
   
    // Exit if invalid parser or EOF.
    if (parser == NULL || parser -> isEOF)
        return;
    
    // If the pointers to nextChromosome and chromosome 
    //  are not equal, then copy the string in nextChromosome
    //  to chromosome.
    if (parser -> nextChrom != chrom) {
        chrom -> l = 0;
        parser -> nextChrom -> l = 0;
        kputs(ks_str(parser -> nextChrom), chrom);
    }
    // Copy all the values from the primed read into the arguments.
    *pos = parser -> nextPos;
    *numOfAlleles = parser -> nextNumAlleles;
    GENOTYPE* temp = *genos;
    *genos = parser -> nextGenos;
    parser -> nextGenos = temp;
    
    // If EOF, set flag and exit.
    if (ks_eof(parser -> stream)) {
        parser -> isEOF = true;
        return;
    }

    // Read the next record into the buffer.
    int dret;
    ks_getuntil(parser -> stream, '\n', parser -> buffer, &dret);

    // If empty line, set EOF and exit.
    if (ks_len(parser -> buffer) == 0) {
        parser -> isEOF = true;
        return;
    }

    // Prime the next read by parsing the record.
    //  Iterate through the buffer string.
    int numTabs = 0, prevIndex = 0, numAlleles = 2;
    for (int i = 0; i <= ks_len(parser -> buffer); i++) {
        if (i == ks_len(parser -> buffer) || ks_str(parser -> buffer)[i] == '\t') {
            // In the first field, set the value of nextChromosome.
            if (numTabs == 0)
                kputsn(ks_str(parser -> buffer), i - prevIndex, parser -> nextChrom);
            // If the second field, set the value of nextPosition.
            else if (numTabs == 1)
                parser -> nextPos = (int) strtol(ks_str(parser -> buffer) + prevIndex + 1, (char**) NULL, 10);
            // If the fifth field, count the number of alleles.
            else if (numTabs == 4) {
                for (int j = prevIndex + 1; ks_str(parser -> buffer)[j] != '\t'; j++)
                    if (ks_str(parser -> buffer)[j] == ',')
                        numAlleles++;
            // If the 9th field or greater, parse each sample's genotype.
            } else if (numTabs > 8)
                parser -> nextGenos[numTabs - 9] = parse_genotype(ks_str(parser -> buffer) + prevIndex + 1, numAlleles);
            prevIndex = i;
            numTabs++;
        }
    }

    parser -> nextNumAlleles = numAlleles;

}

void destroy_vcf_genotype_parser(VCFGenotypeParser* parser) {
    if (parser == NULL)
        return;

    // Free everything sued to read in the file.
    gzclose(parser -> file);
    ks_destroy(parser -> stream);
    free(ks_str(parser -> fileName)); free(parser -> fileName);
    // Free all sample names array.
    for (int i = 0; i < parser -> numSamples; i++)
        free(parser -> sampleNames[i].s);
    free(parser -> sampleNames);
    // Free the buffer.
    free(ks_str(parser -> buffer)); free(parser -> buffer);
    // Free the nextChrom string.
    free(ks_str(parser -> nextChrom)); free(parser -> nextChrom);
    // Free the genotypes array array.
    free(parser -> nextGenos);
    // Free the structure.
    free(parser);
}

// Used to test the parser.
/*
int main() {
    
    VCFGenotypeParser* parser = init_vcf_genotype_parser("vcf_parser_test.vcf.gz");

    printf("There are %d samples with the following names:\n", parser -> num_samples);
    for (int i = 0; i < parser -> numSamples; i++)
        printf("%s\n", ks_str(&(parser -> sampleNames[i])));
    
    kstring_t* chrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    int pos;
    int numOfAlleles;
    GENOTYPE* genos = (GENOTYPE*) calloc(parser -> numSamples, sizeof(GENOTYPE));

    while(true) {

        get_next_locus(parser, chrom, &pos, &numOfAlleles, &genos);

        if (parser -> isEOF)
            break;
        
        printf("\n");
        printf("Chromosome: %s\n", ks_str(chrom));
        printf("Position: %d\n", pos);
        printf("Num Alleles: %d\n", numOfAlleles);
        printf("Genotypes:\n");
        for (int i = 0; i < parser -> numSamples; i++)
            printf("%x\n", genos[i]);
        printf("\n");
        
    }

    destroy_vcf_genotype_parser(parser);
    free(genos);
    free(ks_str(chrom)); free(chrom);
    
}
*/