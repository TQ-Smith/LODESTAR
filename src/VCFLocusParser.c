
// File: 
// Date: 
// Author: TQ Smith
// Purpose: 

#include <stdio.h>

#include "VCFLocusParser.h"

VCFLocusParser* init_vcf_locus_parser(char* fileName, kstring_t* regions, bool takeComplement, double maf, double afMissing, bool dropMonomorphicSites) {

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
                ks_overwriten(ks_str(buffer) + i + 1, i - prevIndex, &sampleNames[numTabs - 8]);
            prevIndex = i;
            numTabs++;
        }

    VCFLocusParser* parser = (VCFLocusParser*) malloc(sizeof(VCFLocusParser));
    parser -> fileName = (kstring_t*) calloc(1, sizeof(kstring_t));
    ks_overwrite(fileName, parser -> fileName);
    parser -> file = file;
    parser -> stream = stream;
    parser -> numSamples = numSamples;
    parser -> sampleNames = sampleNames;
    parser -> buffer = buffer;
    parser -> isEOF = false;
    parser -> nextChrom = (kstring_t*) calloc(1, sizeof(kstring_t));
    parser -> nextLocus = (Locus*) calloc(numSamples, sizeof(Locus));
    if (regions == NULL)
        parser -> filter = NULL;
    else {
        parser -> filter = init_region_filter(regions, takeComplement);
        if (parser -> filter == NULL) {
            destroy_vcf_locus_parser(parser);
            return NULL;
        }
    }
    parser -> maf = maf;
    parser -> afMissing = afMissing;
    parser -> dropMonomorphicSites = dropMonomorphicSites;
    for (int i = 0; i < 16; i++)
        parser -> alleleCounts[i] = 0;

    get_next_locus(parser, parser -> nextChrom, &(parser -> nextPos), &(parser -> nextNumAlleles), &(parser -> nextLocus));

    return parser;

}

void seek(VCFLocusParser* parser) {

    int dret, numTabs, prevIndex, numAlleles;
    bool isInFilter;
    Locus l;
    double maf, afMissing, afMax, af;

    while (true) {

        ks_getuntil(parser -> stream, '\n', parser -> buffer, &dret);

        if (ks_eof(parser -> stream) || ks_len(parser -> buffer) == 0) {
            parser -> isEOF = true;
            return;
        }

        numTabs = 0, prevIndex = 0, numAlleles = 2;
        isInFilter = true;
        for (int i = 0; i <= ks_len(parser -> buffer) && isInFilter; i++) {
            if (i == ks_len(parser -> buffer) || ks_str(parser -> buffer)[i] == '\t') {
                if (numTabs == 0) {
                    ks_overwriten(ks_str(parser -> buffer), i, parser -> nextChrom);
                } else if (numTabs == 1) {
                    parser -> nextPos = (int) strtol(ks_str(parser -> buffer) + prevIndex + 1, (char**) NULL, 10);
                    if (parser -> filter != NULL && !query_locus(parser -> filter, parser -> nextChrom, (unsigned int) parser -> nextPos))
                        isInFilter = false;
                } else if (numTabs == 4) {
                    for (int j = prevIndex + 1; ks_str(parser -> buffer)[j] != '\t'; j++)
                        if (ks_str(parser -> buffer)[j] == ',')
                            numAlleles++;
                } else if (numTabs > 8) {
                    l = parse_locus(ks_str(parser -> buffer) + prevIndex + 1, numAlleles);
                    parser -> alleleCounts[LEFT_ALLELE(l)]++;
                    parser -> alleleCounts[RIGHT_ALLELE(l)]++;
                    parser -> nextLocus[numTabs - 9] = l;
                }
                prevIndex = i;
                numTabs++;
            }
        }

        if (!isInFilter)
            continue;
        
        afMissing = parser -> alleleCounts[numAlleles] / (2.0 * parser -> numSamples);
        parser -> alleleCounts[numAlleles] = 0;
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
        
        parser -> nextNumAlleles = numAlleles;

        if (parser -> dropMonomorphicSites && afMax == 1)
            continue;

        if (maf >= parser -> maf && afMissing < parser -> afMissing)
            return;

    }

}

void get_next_locus(VCFLocusParser* parser, kstring_t* chrom, unsigned int* pos, int* numOfAlleles, Locus** locus) {
    if (parser == NULL || parser -> isEOF)
        return;
    
    if (parser -> nextChrom != chrom) {
        chrom -> l = 0;
        parser -> nextChrom -> l = 0;
        ks_overwrite(ks_str(parser -> nextChrom), chrom);
    }

    *pos = parser -> nextPos;
    *numOfAlleles = parser -> nextNumAlleles;
    Locus* temp = *locus;
    *locus = parser -> nextLocus;
    parser -> nextLocus = temp;

    seek(parser);
}


void destroy_vcf_locus_parser(VCFLocusParser* parser) {
    if (parser == NULL)
        return;
    gzclose(parser -> file);
    ks_destroy(parser -> stream);
    free(ks_str(parser -> fileName)); free(parser -> fileName);
    if (parser -> filter != NULL)
        destroy_region_filter(parser -> filter);
    for (int i = 0; i < parser -> numSamples; i++)
        free(parser -> sampleNames[i].s);
    free(parser -> sampleNames);
    free(ks_str(parser -> buffer)); free(parser -> buffer);
    free(ks_str(parser -> nextChrom)); free(parser -> nextChrom);
    free(parser -> nextLocus);
    free(parser);
}

// Used to test the parser.
/*
int main() {
    
    kstring_t* intervals = (kstring_t*) calloc(1, sizeof(kstring_t));
    ks_overwrite("3", intervals);
    VCFLocusParser* parser = init_vcf_locus_parser("./data/vcf_parser_test.vcf.gz", intervals, false, 0, 1, false);

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