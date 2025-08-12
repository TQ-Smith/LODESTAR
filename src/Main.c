
// File: Main.c
// Date: 
// Version 1: 
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Run LODESTAR analysis.

#include "Interface.h"
#include "LODESTAR.h"
#include <stdio.h>
#include <stdlib.h>

void print_json_matrix(FILE* out, double** X, int N, int K) {
    fprintf(out, "[\n");
    for (int i = 0; i < N; i++) {
        fprintf(out, "[ ");
        for (int j = 0; j < K - 1; j++) {
            fprintf(out, "%lf, ", X[i][j]);
        }
        fprintf(out, "%lf ]", X[i][K - 1]);
        if (i != N - 1)
            fprintf(out, ",\n");
        else
            fprintf(out, "\n");
    }
    fprintf(out, "]");
}

void print_json(LodestarConfig_t* lodestarConfig, BlockList_t* globalList, double** y, double* y0) {
    kstring_t* outName = calloc(1, sizeof(kstring_t));
    ksprintf(outName, "%s.json", lodestarConfig -> outputBasename);
    FILE* out = fopen(outName -> s, "w");
    free(outName -> s); free(outName);

    fprintf(out, "{\n");

    // Print out all the global information.
    fprintf(out, "\"Command\": \"%s\",\n", lodestarConfig -> cmd);
    fprintf(out, "\"GlobalNumberOfLoci\": %d,\n", globalList -> numLoci);
    fprintf(out, "\"GlobalNumberOfHaplotypes\": %d,\n", globalList -> numHaps);
    fprintf(out, "\"GlobalVarainceCaptured\": %lf,\n", globalList -> varCapt);
    if (lodestarConfig -> targetFileName == NULL) {
        fprintf(out, "\"GlobalProcrustesStatistic\": null,\n");
        fprintf(out, "\"GlobalProcrustesStatisticPvalue\": null,\n");
    } else {
        fprintf(out, "\"GlobalProcrustesStatistic\": %lf,\n", globalList -> procrustesT);
        fprintf(out, "\"GlobalProcrustesStatisticPvalue\": %lf,\n", globalList -> pvalue);
    }
    fprintf(out, "\"GlobalX\":");
    print_json_matrix(out, globalList -> X, globalList -> numSamples, lodestarConfig -> k);
    fprintf(out, ",\n");
    if (lodestarConfig -> targetFileName == NULL) {
        fprintf(out, "\"Y\": null,\n");
        fprintf(out, "\"y0\": null,\n");
    } else {
        fprintf(out, "\"Y\":");
        print_json_matrix(out, y, globalList -> numSamples, lodestarConfig -> k);
        fprintf(out, ",\n");
        fprintf(out, "\"y0\": [ ");
        for (int i = 0; i < lodestarConfig -> k - 1; i++)
            fprintf(out, "%lf, ", y0[i]);
        fprintf(out, "%lf ],\n", y0[lodestarConfig -> k - 1]);
    }
    
    // Print out all the blocks.
    fprintf(out, "\t\"Blocks\": [");
    for (Block_t* temp = globalList -> head; temp != NULL; temp = temp -> next) {
        if (temp -> isDropped)
            continue;
        fprintf(out, "\t{\n");
        fprintf(out, "\t\"BlockNumber\": %d,\n", temp -> blockNum);
        fprintf(out, "\t\"BlockNumberOnChromosome\": %d,\n", temp -> blockNumOnChrom);
        fprintf(out, "\t\"Chromosome\": \"%s\",\n", temp -> chrom);
        fprintf(out, "\t\"StartCoordinate\": %d,\n", temp -> startCoordinate);
        fprintf(out, "\t\"EndCoordinate\": %d,\n", temp -> endCoordinate);
        fprintf(out, "\t\"NumberOfLoci\": %d,\n", temp -> numLoci);
        fprintf(out, "\t\"NumberOfHaplotypes\": %d,\n", temp -> numHaps);
        if (temp -> isDropped) {
            fprintf(out, "\t\"VarianceCaptured\": %d,\n", -1);
            fprintf(out, "\t\"ProcrustesStatistic\": %d,\n", -1);
            fprintf(out, "\t\"ProcrustesPValue\": %d,\n", -1);
        } else {
            fprintf(out, "\t\"VarianceCaptured\": %lf,\n", temp -> varCapt);
            fprintf(out, "\t\"ProcrustesStatistic\": %lf,\n", temp -> procrustesT);
            fprintf(out, "\t\"ProcrustesPValue\": %lf,\n", temp -> pvalue);
        }
        fprintf(out, "\t\"X\": ");
        if (temp -> X != NULL)
            print_json_matrix(out, temp -> X, globalList -> numSamples, lodestarConfig -> k);
        else 
            fprintf(out, "NULL\n");
        fprintf(out, "}");
        if (temp -> next != NULL)
            fprintf(out, ",\n");
        else
            fprintf(out, "\n");
    }
    fprintf(out, "\t]\n");
    fprintf(out, "}\n");

    fclose(out);
}

void print_summary(LodestarConfig_t* lodestarConfig, BlockList_t* globalList) {
    kstring_t* outName = calloc(1, sizeof(kstring_t));
    ksprintf(outName, "%s.tsv", lodestarConfig -> outputBasename);
    FILE* out = fopen(outName -> s, "w");
    free(outName -> s); free(outName);

    // Echo command.
    fprintf(out, "#%s\n", lodestarConfig -> cmd);

    // Print header.
    fprintf(out, "BlockNum\tBlockNumOnChr\tChr\tStart\tEnd\tNumLoci\tNumHaps\tVarianceCaptured\tProcrustesStatistic\tP-Value\n");

    // Print out each block.
    for(Block_t* temp = globalList -> head; temp != NULL; temp = temp -> next) {
        if (temp -> isDropped)
            continue;
        fprintf(out, "%d\t%d\t%s\t%d\t%d\t%d\t%d\t", temp -> blockNum, temp -> blockNumOnChrom, temp -> chrom, temp -> startCoordinate, temp -> endCoordinate, temp -> numLoci, temp -> numHaps);
        fprintf(out, "%lf\t%lf\t%lf\n", temp -> varCapt, temp -> procrustesT, temp -> pvalue);
    }
    fprintf(out, "0\t0\tGLOBAL\t0\t0\t%d\t%d\t%lf\t", globalList -> numLoci, globalList -> numHaps, globalList -> varCapt);
    if (globalList -> procrustesT != -1)
        fprintf(out, "%lf\t%lf\n", globalList -> procrustesT, globalList -> pvalue);
    else 
        fprintf(out, "-1\t-1\n");

    fclose(out);
}

double** open_target_file(char* targetFileName, int N, int K) {
    // Create target matrix.
    double** targetPoints = init_matrix(N, K);

    // Read in values from target.
    int numLines = 0;
    size_t length = 0;
    char* line = calloc(1024, sizeof(char));
    FILE* targetFile = fopen(targetFileName, "r");
    while (getline(&line, &length, targetFile) != -1) {
        char* tok = strtok(line, "\t");
        for (int i = 0; tok != NULL && i < K; i++) {
            targetPoints[numLines][i] = atof(tok);
            tok = strtok(NULL, "\t");
        }
        // If invalid formatting, free all memory and exit.
        if (tok != NULL || numLines == N) {
            fclose(targetFile);
            free(line);
            destroy_matrix(targetPoints, N);
            return NULL;
        }
        numLines++;
    }
    fclose(targetFile);
    free(line);
    return targetPoints;
}


int main (int argc, char *argv[]) {

    // Get configuration.
    LodestarConfig_t* lodestarConfig = init_lodestar_config(argc, argv);
    if (lodestarConfig == NULL)
        return -1;
   
    // We create the VCF parser and haplotype encoder.
    VCFLocusParser_t* parser = init_vcf_locus_parser(lodestarConfig -> inputFileName, lodestarConfig -> maf, lodestarConfig -> afMissing, true);
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(parser -> numSamples);
    
    // Check to make sure k is valid.
    if (lodestarConfig -> k >= parser -> numSamples) {
        fprintf(stderr, "k = %d >= numSamples %d. Exiting!\n", lodestarConfig -> k, parser -> numSamples);
        free(lodestarConfig);
        destroy_vcf_locus_parser(parser);
        return -1;
    }

    // User defined points to perform Procrustes against.
    double** y = NULL;
    double* y0 = NULL;

    // If target file supplied, make sure it is valid, at least numSamples-by-k
    if (lodestarConfig -> targetFileName != NULL) {
        y = open_target_file(lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
        if (y == NULL) {
            fprintf(stderr, "-y %s is not an %d-by-%d matrix. Exiting!\n", lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
            free(lodestarConfig);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
        // Center and normalize.
        y0 = calloc(lodestarConfig -> k, sizeof(double));
        center_matrix(y, y0, parser -> numSamples, lodestarConfig -> k);
        normalize_matrix(y, parser -> numSamples, lodestarConfig -> k);
    }
    
    // Partition genome into blocks and calculate IBS within the blocks.
    fprintf(stderr, "Beginning blocking algorithm ...\n");
    BlockList_t* globalList = block_allele_sharing(parser, encoder, encoder -> numSamples, lodestarConfig -> BLOCK_SIZE, lodestarConfig -> HAP_SIZE, lodestarConfig -> dropThreshold, lodestarConfig -> threads);

    // Convert IBS to ASD and compute MDS on global.
    fprintf(stderr, "\nFinished blocking algorithm. Starting genome-wide MDS calculations ...\n");
    double* asd = calloc(PACKED_SIZE(encoder -> numSamples), sizeof(double));
    RealSymEigen_t* eigen = init_real_sym_eigen(encoder -> numSamples);
    for (int i = 0; i < encoder -> numSamples; i++)
        for (int j = i + 1; j < encoder -> numSamples; j++)
            asd[PACKED_INDEX(i, j)] = ibs_to_asd(globalList -> alleleCounts[PACKED_INDEX(i, j)]);
    globalList -> X = init_matrix(encoder -> numSamples, lodestarConfig -> k);
    globalList -> varCapt = compute_classical_mds(eigen, asd, lodestarConfig -> k, globalList -> X);
    if (y != NULL)
        globalList -> procrustesT = procrustes_statistic(globalList -> X, NULL, y, y0, eigen, eigen -> N, lodestarConfig -> k, true);
    destroy_real_sym_eigen(eigen);
    free(asd);

    // Convert IBS to ASD and calculate jackknifed procrustes statistic.
    fprintf(stderr, "Finished genome-wide MDS calulations. Starting Procrustes ...\n\n");
    if (lodestarConfig -> sampleSize == 0)
        lodestarConfig -> sampleSize = globalList -> numBlocks;
    procrustes(globalList, y, y0, lodestarConfig -> k, lodestarConfig -> threads, lodestarConfig -> numReps, lodestarConfig -> sampleSize);
    if (lodestarConfig -> numReps == 0)
        fprintf(stderr, "Writing results to output files...\n");
    else
        fprintf(stderr, "\nFinished Bootstrap. Writing results to output files...\n");

    // Print summary and JSON file.
    print_summary(lodestarConfig, globalList);
    print_json(lodestarConfig, globalList, y, y0);

    // Free used memory.
    if (y0 != NULL) free(y0);
    destroy_matrix(y, encoder -> numSamples);
    destroy_block_list(globalList);
    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);
    destroy_lodestar_config(lodestarConfig);
    fprintf(stderr, "Done!\n");
    return 0;
}