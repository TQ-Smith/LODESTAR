
// File: 
// Date: 
// Author: TQ Smith
// Purpose: 

#include <stdio.h>

#include "VCFLocusParser.h"

VCFLocusParser* init_vcf_locus_parser(char* fileName) {

    gzFile file = gzopen(fileName, "r");

    int errnum;
    gzerror(file, &errnum);
    if (errnum != Z_OK)
        return NULL;
    
    kstream_t* stream = ks_init(file);
    
    kstring_t* buffer = (kstring_t*) calloc(1, sizeof(kstring_t));
    int dret;

    do {
        ks_getuntil(stream, '\n', buffer, &dret);
    } while (strncmp(ks_str(buffer), "#C", 2) != 0);
    
    int numSamples = 0;
    for (int i = 0; i < ks_len(buffer); i++)
        if (buffer -> s[i] == '\t')
            numSamples++;
    numSamples -= 8;

    kstring_t* sampleNames = (kstring_t*) calloc(numSamples, sizeof(kstring_t));
    int numTabs = 0, prevIndex;
    for (int i = 0; i <= ks_len(buffer); i++)
        if (ks_str(buffer)[i] == '\t') {
            if (numTabs > 7)
                kputsn(ks_str(buffer) + i + 1, i - prevIndex, &sampleNames[numTabs - 8]);
            prevIndex = i;
            numTabs++;
        }

    VCFLocusParser* parser = (VCFLocusParser*) malloc(sizeof(VCFLocusParser));
    parser -> fileName = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(fileName, parser -> fileName);
    parser -> file = file;
    parser -> stream = stream;
    parser -> numSamples = numSamples;
    parser -> sampleNames = sampleNames;
    parser -> buffer = buffer;
    parser -> isEOF = false;
    parser -> nextChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    parser -> nextLoci = (Locus*) calloc(numSamples, sizeof(Locus));

    get_next_locus(parser, parser -> nextChrom, &(parser -> nextPos), &(parser -> nextNumAlleles), &(parser -> nextLoci));

    return parser;

}

void get_next_locus(VCFLocusParser* parser, kstring_t* chrom, int* pos, int* numOfAlleles, Locus** loci) {
   
    if (parser == NULL || parser -> isEOF)
        return;
    
    if (parser -> nextChrom != chrom) {
        chrom -> l = 0;
        parser -> nextChrom -> l = 0;
        kputs(ks_str(parser -> nextChrom), chrom);
    }

    *pos = parser -> nextPos;
    *numOfAlleles = parser -> nextNumAlleles;
    Locus* temp = *loci;
    *loci = parser -> nextLoci;
    parser -> nextLoci = temp;
    
    if (ks_eof(parser -> stream)) {
        parser -> isEOF = true;
        return;
    }

    int dret;
    ks_getuntil(parser -> stream, '\n', parser -> buffer, &dret);

    if (ks_len(parser -> buffer) == 0) {
        parser -> isEOF = true;
        return;
    }

    int numTabs = 0, prevIndex = 0, numAlleles = 2;
    for (int i = 0; i <= ks_len(parser -> buffer); i++) {
        if (i == ks_len(parser -> buffer) || ks_str(parser -> buffer)[i] == '\t') {
            if (numTabs == 0)
                kputsn(ks_str(parser -> buffer), i - prevIndex, parser -> nextChrom);
            else if (numTabs == 1)
                parser -> nextPos = (int) strtol(ks_str(parser -> buffer) + prevIndex + 1, (char**) NULL, 10);
            else if (numTabs == 4) {
                for (int j = prevIndex + 1; ks_str(parser -> buffer)[j] != '\t'; j++)
                    if (ks_str(parser -> buffer)[j] == ',')
                        numAlleles++;
            } else if (numTabs > 8)
                parser -> nextLoci[numTabs - 9] = parse_locus(ks_str(parser -> buffer) + prevIndex + 1, numAlleles);
            prevIndex = i;
            numTabs++;
        }
    }

    parser -> nextNumAlleles = numAlleles;

}

void destroy_vcf_locus_parser(VCFLocusParser* parser) {
    if (parser == NULL)
        return;

    gzclose(parser -> file);
    ks_destroy(parser -> stream);
    free(ks_str(parser -> fileName)); free(parser -> fileName);
    for (int i = 0; i < parser -> numSamples; i++)
        free(parser -> sampleNames[i].s);
    free(parser -> sampleNames);
    free(ks_str(parser -> buffer)); free(parser -> buffer);
    free(ks_str(parser -> nextChrom)); free(parser -> nextChrom);
    free(parser -> nextLoci);
    free(parser);
}
