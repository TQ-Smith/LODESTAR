
// File: main.c
// Date: 
// Version 1: 
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Command-line interface for LODESTAR and main analysis.

#include "Interface.h"
#include "SlidingWindow.h"
#include "ProcrustesAnalysis.h"
#include "Logger.h"
#include "Matrix.h"
MATRIX_INIT(double, double)
#include <time.h>

int main (int argc, char *argv[]) {

    // Seed random number generator.
    srand(time(NULL));

    // Get configuration.
    LodestarConfiguration_t* lodestar_config = init_lodestar_config(argc, argv);
    if (lodestar_config == NULL) {
        fprintf(stderr, "Exiting!\n");
        return -1;
    }

    // We create the VCF parser.
    VCFLocusParser_t* parser = init_vcf_locus_parser(lodestar_config -> inputFileName, false, lodestar_config -> maf, lodestar_config -> afMissing, true);

    // Check that VCF file is valid.
    if (parser == NULL) {
        fprintf(stderr, "Invalid VCF file. Exiting.\n");
        free(lodestar_config);
        return -1;
    }

    // Check to make sure k is valid.
    if (lodestar_config -> k >= parser -> numSamples) {
        fprintf(stderr, "k = %d >= numSamples %d. Exiting.\n", lodestar_config -> k, parser -> numSamples);
        free(lodestar_config);
        destroy_vcf_locus_parser(parser);
        return -1;
    }

    // User defined points to perform Procrustes against.
    double** userPoints = NULL;

    // If target file supplied, make sure it is valid, at least numSamples-by-k
    if (lodestar_config -> targetFileName != NULL) {
        userPoints = open_target_file(lodestar_config -> targetFileName, parser -> numSamples, lodestar_config -> k);
        if (userPoints == NULL) {
            fprintf(stderr, "--target %s does not contain %d values. Exiting.\n", lodestar_config -> targetFileName, lodestar_config -> k * (parser -> numSamples));
            free(lodestar_config);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
    }

    // Our group information.
    int* samplesToGroups = NULL;
    int numGroups = 0;

    // If group file supplied, make sure it is valid. numSamples-by-2
    if (lodestar_config -> groupFileName != NULL) {
        samplesToGroups = open_group_file(lodestar_config -> groupFileName, parser -> sampleNames, parser -> numSamples, &numGroups);
        if (samplesToGroups == NULL) {
            fprintf(stderr, "--group %s does not contain sample-group pairs. Exiting.\n", lodestar_config -> targetFileName);
            free(lodestar_config);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
    }

    // Create the haplotype encoder.
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(parser -> numSamples);

    // Valid configuration and input. We being our LODESTAR analysis.

    // Create output files.
    FILE* windowSummaries = NULL;
    FILE* windowPoints = NULL;
    FILE* windowTarget = NULL;
    kstring_t* outputBasename = init_kstring(lodestar_config -> outputBasename);
    printf("\nLogging progress in %s.log\n\n", ks_str(outputBasename));
    kputs(".log", outputBasename);
    INIT_LOG(ks_str(outputBasename));
    ks_overwrite(lodestar_config -> outputBasename, outputBasename);
    if (!lodestar_config -> global) {
        kputs("_summary.tsv", outputBasename);
        windowSummaries = fopen(ks_str(outputBasename), "w");
        ks_overwrite(lodestar_config -> outputBasename, outputBasename);
    }
    // If we are using JSON output.
    if (lodestar_config -> useJsonOutput) {
        kputs("_windows.json", outputBasename);
        windowPoints = fopen(ks_str(outputBasename), "w");
    } else {
        kputs("_windows.txt", outputBasename);
        windowPoints = fopen(ks_str(outputBasename), "w");
    }
    ks_overwrite(lodestar_config -> outputBasename, outputBasename);
    kputs("_target.tsv", outputBasename);
    windowTarget = fopen(ks_str(outputBasename), "w");
    destroy_kstring(outputBasename);

    // Print initial configuration to log file.
    print_configuration(logger -> file, lodestar_config);
    fprintf(logger -> file, "\nSample Names:\n");
    for (int i = 0; i < parser -> numSamples; i++)
        fprintf(logger -> file, "%s\n", ks_str(parser -> sampleNames[i]));
    fprintf(logger -> file, "\n");

    Window_t* global = NULL;
    Window_t** windows = NULL;
    int numWindows = 0;

    if (lodestar_config -> global) {
        LOG_INFO("Performing genome-wide calculations ...\n");
        printf("Performing genome-wide calculations ...\n\n");
        global = global_window(parser, encoder, lodestar_config -> k, lodestar_config -> HAP_SIZE, lodestar_config -> threads);
        LOG_INFO("Finished genome-wide calculations ...\n");
        printf("Finished genome-wide calculations ...\n\n");
    } else {
        LOG_INFO("Performing sliding-window calculations ...\n");
        printf("Performing sliding-window calculations ...\n\n");
        windows = sliding_window(parser, encoder, lodestar_config -> k, lodestar_config -> HAP_SIZE, lodestar_config -> STEP_SIZE, lodestar_config -> WINDOW_SIZE, lodestar_config -> threads, lodestar_config -> MAX_GAP, &numWindows);
        LOG_INFO("Finished sliding-window calculations ...\n");
        printf("Finished sliding-window calculations ...\n\n");
    }

    // We now configure the target to perform Procrustes analysis.
    double** targetPoints = NULL;
    double* targetPointsColMeans = NULL;

    // Group by user coordinates.
    if (samplesToGroups != NULL && userPoints != NULL) {
        targetPoints = assign_centroid_of_group(userPoints, samplesToGroups, parser -> numSamples, lodestar_config -> k, numGroups);
        targetPointsColMeans = (double*) calloc(lodestar_config -> k, sizeof(double));
        center_matrix(targetPoints, targetPointsColMeans, parser -> numSamples, lodestar_config -> k);
        // Destroy our user points. Not needed anymore.
        destroy_matrix(double, userPoints, parser -> numSamples);
    // Group by global coordinates.
    } else if (samplesToGroups != NULL && windows != NULL && windows[0] -> X != NULL) {
        targetPoints = assign_centroid_of_group(windows[0] -> X, samplesToGroups, parser -> numSamples, lodestar_config -> k, numGroups);
        targetPointsColMeans = (double*) calloc(lodestar_config -> k, sizeof(double));
        center_matrix(targetPoints, targetPointsColMeans, parser -> numSamples, lodestar_config -> k);
    // Use the user supplied coordinates.
    } else if (userPoints != NULL) {
        targetPoints = userPoints;
        // Mean center input target points.
        targetPointsColMeans = (double*) calloc(lodestar_config -> k, sizeof(double));
        center_matrix(targetPoints, targetPointsColMeans, parser -> numSamples, lodestar_config -> k);
    // Use the global coordinates.
    } else if (windows[0] -> X != NULL) {
        targetPoints = windows[0] -> X;
        targetPointsColMeans = NULL;
    // Otherwise, error. We do not have a valid matrix to perform Procrustes against.
    } else if (!lodestar_config -> global) {
        LOG_WARNING("Unable to perform Procrustes. Either no user supplied points were given or global did not converge!\n");
    }

    // Used for Procrustes analysis and to transform points.
    RealSymEigen_t* eigen = init_real_sym_eigen(parser -> numSamples);

    // We only perform Procrustes if the target points was set.
    if (targetPoints != NULL) {
        LOG_INFO("Beginning Procrustes Analysis ...\n");
        printf("Beginning Procrustes Analysis ...\n\n");

        // Just global against the target.
        if (global != NULL) {
            double** shuffleX = create_matrix(double, parser -> numSamples, lodestar_config -> k);
            global -> t = procrustes_statistic(global -> X, NULL, targetPoints, NULL, eigen, eigen -> N, lodestar_config -> k, false, lodestar_config -> similarity);
            global -> pval = permutation_test(global -> X, targetPoints, shuffleX, eigen, eigen -> N, lodestar_config -> k, lodestar_config -> similarity, global -> t, lodestar_config -> NUM_PERMS);
            // Trasnsform global set of points.
            global -> t = procrustes_statistic(global -> X, NULL, targetPoints, targetPointsColMeans, eigen, eigen -> N, lodestar_config -> k, true, lodestar_config -> similarity);
            destroy_matrix(double, shuffleX, encoder -> numSamples);
        // Sliding window against the target.
        } else {
            procrustes_sliding_window(windows, numWindows, targetPoints, parser -> numSamples, lodestar_config -> k, lodestar_config -> similarity, lodestar_config -> pthresh == 0 ? 0 : lodestar_config -> NUM_PERMS, lodestar_config -> threads);
        }

        LOG_INFO("Finished Procrustes Analysis ...\n");
        printf("Finished Procrustes Analysis ...\n\n");
    }

    LOG_INFO("Saving results to output files ...\n");
    printf("Saving results to output files ...\n\n");

    // Output results.
    if (targetPoints == NULL) {
        if (lodestar_config -> useJsonOutput) fprintf(windowPoints, "{\n\t\"windows\": [");
        print_window(windowPoints, parser -> sampleNames, windows[0], parser -> numSamples, lodestar_config -> k, lodestar_config -> useJsonOutput, true);
        if (lodestar_config -> useJsonOutput) fprintf(windowPoints, "\n\t]\n}\n");
    } else {
        // Echo command in summary file for convience.
        fprintf(windowSummaries, "#Command: ");
        for (int i = 0; i < argc; i++) 
            fprintf(windowSummaries, "%s ", argv[i]);
        fprintf(windowSummaries, "\nWin\tWinChr\tChr\tStart\tEnd\tnLoci\tnHaps\tp-val\tt-stat\n");
        if (lodestar_config -> useJsonOutput) fprintf(windowPoints, "{\n\t\"windows\": [");
        for (int i = 1; i < numWindows; i++) {
            print_window_summary(windowSummaries, windows[i]);
            fprintf(windowPoints, "\n");
            // Print points if thersholds are met. Print window points.
            if ((windows[i] -> X != NULL) && ((lodestar_config -> pthresh != 0 && windows[i] -> pval < lodestar_config -> pthresh) || (windows[i] -> t >= lodestar_config -> tthresh))) {
                procrustes_statistic(windows[i] -> X, NULL, targetPoints, targetPointsColMeans, eigen, parser -> numSamples, lodestar_config -> k, true, lodestar_config -> similarity);
                print_window(windowPoints, parser -> sampleNames, windows[i], parser -> numSamples, lodestar_config -> k, lodestar_config -> useJsonOutput, true);
            } else {
                print_window(windowPoints, parser -> sampleNames, windows[i], parser -> numSamples, lodestar_config -> k, lodestar_config -> useJsonOutput, false);
            }
            if (lodestar_config -> useJsonOutput) fprintf(windowPoints, ",");
        }
        // Print the global window.
        fprintf(windowPoints, "\n");
        print_window(windowPoints, parser -> sampleNames, windows[0], parser -> numSamples, lodestar_config -> k, lodestar_config -> useJsonOutput, true);
        if (lodestar_config -> useJsonOutput) fprintf(windowPoints, "\n\t]\n}\n");

        // Print the target used.
        for (int i = 0; i < parser -> numSamples; i++) {
            for (int j = 0; j < lodestar_config -> k - 1; j++) {
                fprintf(windowTarget, "%lf\t", targetPoints[i][j]);
            }
            fprintf(windowTarget, "%lf\n", targetPoints[i][lodestar_config -> k - 1]);
        }
    }

    LOG_INFO("Finished Analysis! Exiting ...\n");
    printf("Finished Analysis! Exiting ...\n\n");

    // Close all files.
    if (windowSummaries)
        fclose(windowSummaries);
    if (windowPoints)
        fclose(windowPoints);
    if (windowTarget)
        fclose(windowTarget);
    CLOSE_LOG();

    // Free all used memory.
    if (targetPoints != NULL && targetPoints != windows[0] -> X)
        destroy_matrix(double, targetPoints, parser -> numSamples);
    if (targetPointsColMeans != NULL)
        free(targetPointsColMeans);
    if (windows != NULL) {
        for (int i = 0; i < numWindows; i++)
            destroy_window(windows[i], parser -> numSamples);
        free(windows);
    }
    if (global != NULL) {
        destroy_window(global, parser -> numSamples);
    }
    if (samplesToGroups)
        free(samplesToGroups);
    free(lodestar_config);
    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);
    destroy_real_sym_eigen(eigen);
    return 0;
}