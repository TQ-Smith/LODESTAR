
// File: Interface.c
// Date: 18 February 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines LODESTAR configuration from command line interface and
//              output formatting of LODESTAR results.

#include "Interface.h"
#include <unistd.h> // Check if file exits.
#include "../lib/ketopt.h"
#include "Matrix.h"
MATRIX_INIT(double, double)

double** open_target_file(char* targetFileName, int N, int K) {
    double inVal;
    int numVals = 0;
    // Create target matrix.
    double** targetPoints = create_matrix(double, N, K);
    // Read in values from target.
    FILE* targetFile = fopen(targetFileName, "r");
    while (numVals < N * K && fscanf(targetFile, "%lf", &inVal) != EOF) {
        targetPoints[numVals / K][numVals % K] = inVal;
        numVals++;
    }
    fclose(targetFile);
    // If file is invalid, then print error and free used memory.
    if (numVals != N * K) {
        destroy_matrix(double, targetPoints, N);
        return NULL;
    }
    return targetPoints;
}

#define PRINT_BOOL(X) (X ? "true" : "false")

// Print verbose LODESTAR configuration.
// Accepts:
//  FILE* output -> The output stream to print the configuration.
//  LodestarConfiguration_t* lodestart_config -> The configuration to print.
// Returns: void.
void print_configuration(FILE* output, LodestarConfiguration_t* lodestar_config) {
    fprintf(output, "LODESTAR Configuration\n");
    fprintf(output, "----------------------\n");
    fprintf(output, "Input file: %s\n", lodestar_config -> inputFileName);
    fprintf(output, "Output basename: %s\n", lodestar_config -> outputBasename);
    fprintf(output, "Calculate genome-wide only: %s\n", PRINT_BOOL(lodestar_config -> global));
    fprintf(output, "Haplotype size: %d\n", lodestar_config -> HAP_SIZE);
    fprintf(output, "Window size: %d\n", lodestar_config -> WINDOW_SIZE);
    fprintf(output, "Step size: %d\n", lodestar_config -> STEP_SIZE);
    fprintf(output, "K: %d\n", lodestar_config -> k);
    fprintf(output, "Use similarity: %s\n", PRINT_BOOL(lodestar_config -> similarity));
    fprintf(output, "Number of threads used: %d\n", lodestar_config -> threads);
    fprintf(output, "P-value threshold (Disabled if -1): %lf\n", lodestar_config -> pthresh);
    fprintf(output, "Number of permutations (Disabled when p-value is -1): %d\n", lodestar_config -> NUM_PERMS);
    fprintf(output, "Procrustes statistic threshold: %lf\n", lodestar_config -> tthresh);
    fprintf(output, "Biallelic minor allele frequency threshold: %lf\n", lodestar_config -> maf);
    fprintf(output, "Missing allele frequency threshold: %lf\n", lodestar_config -> afMissing);
    fprintf(output, "Max range of window in basepairs: %d\n", lodestar_config -> MAX_GAP);
    if (lodestar_config -> targetFileName != NULL)
        fprintf(output, "File of coordinates to perform Procrustes analysis against: %s\n", lodestar_config -> targetFileName);
    fprintf(output, "Save as JSON file: %s\n", PRINT_BOOL(lodestar_config -> useJsonOutput));
}

void print_window_summary(FILE* output, Window_t* window) {
    fprintf(output, "%d\t", window -> winNum);
    fprintf(output, "%d\t", window -> winNumOnChrom);
    fprintf(output, "%s\t", ks_str(window -> chromosome));
    fprintf(output, "%d\t", window -> startCoord);
    fprintf(output, "%d\t", window -> endCoord);
    fprintf(output, "%d\t", window -> numLoci);
    fprintf(output, "%d\t", window -> numHaps);
    if (window -> t == -1) {
        fprintf(output, "-1\t");
        fprintf(output, "-1\n");
    } else {
        fprintf(output, "%lf\t", window -> pval);
        fprintf(output, "%lf\n", window -> t);
    }
}

// Print a row of a matrix.
// Accepts:
//  FILE* output -> The output stream to print to.
//  double* row -> The array to print.
//  int K -> The number of elements of row to print.
//  bool useJsonOutput -> Flag to indicate if row should be JSON formatted.
// Returns: void.
void print_row(FILE* output, double* row, int K, bool useJsonOutput) {
    if (useJsonOutput) {
        for (int i = 0; i < K - 1; i++)
            fprintf(output, "%lf, ", row[i]);
        fprintf(output, "%lf]", row[K - 1]);
    } else {
        for (int i = 0; i < K - 1; i++)
            fprintf(output, "%lf\t", row[i]);
        fprintf(output, "%lf", row[K - 1]);
    }
}

void print_window(FILE* output, kstring_t** sampleNames, Window_t* window, int N, int K, bool useJsonOutput, bool printCoords) {
    // Print JSON.
    if (useJsonOutput) {
        // Print single fields.
        fprintf(output, "\t\t{\n");
        fprintf(output, "\t\t\t\"Window Number\": %d,\n", window -> winNum);
        fprintf(output, "\t\t\t\"Window Number on Chromosome\": %d,\n", window -> winNumOnChrom);
        fprintf(output, "\t\t\t\"Chromosome\": \"%s\",\n", ks_str(window -> chromosome));
        fprintf(output, "\t\t\t\"Start Coordinate\": %d,\n", window -> startCoord);
        fprintf(output, "\t\t\t\"End Coordinate\": %d,\n", window -> endCoord);
        fprintf(output, "\t\t\t\"Number of Loci\": %d,\n", window -> numLoci);
        fprintf(output, "\t\t\t\"Number of Haplotypes\": %d,\n", window -> numHaps);
        if (window -> t == -1) {
            fprintf(output, "\t\t\t\"P-Value\": -1,\n");
            fprintf(output, "\t\t\t\"t-statistic\": -1,\n");
        } else {
            fprintf(output, "\t\t\t\"P-Value\": %lf,\n", window -> pval);
            fprintf(output, "\t\t\t\"t-statistic\": %lf,\n", window -> t);
        }
        // Print points.
        fprintf(output, "\t\t\t\"Points\": ");
        if (!printCoords || window -> X == NULL) {
            fprintf(output, "null\n");
        } else {
            fprintf(output, "[\n");
            for (int i = 0; i < N; i++) {
                fprintf(output, "\t\t\t\t[");
                print_row(output, window -> X[i], K, true);
                if (i != N - 1)
                    fprintf(output, ",\n");
                else
                    fprintf(output, "\n");
            }
            fprintf(output, "\t\t\t]\n");
        }
        fprintf(output, "\t\t}");
    // Print as plain text.
    } else {
        // Print single fields.
        fprintf(output, "Window Number: %d\n", window -> winNum);
        fprintf(output, "Window Number on Chromosome: %d\n", window -> winNumOnChrom);
        fprintf(output, "Chromosome: %s\n", ks_str(window -> chromosome));
        fprintf(output, "Start Coordinate: %d\n", window -> startCoord);
        fprintf(output, "End Coordinate: %d\n", window -> endCoord);
        fprintf(output, "Number of Loci: %d\n", window -> numLoci);
        fprintf(output, "Number of Haplotypes: %d\n", window -> numHaps);
        if (window -> t == -1) {
            fprintf(output, "P-Value: -1\n");
            fprintf(output, "t-statistic: -1\n");
        } else {
            fprintf(output, "P-Value: %lf\n", window -> pval);
            fprintf(output, "t-statistic: %lf\n", window -> t);
        }
        // Print points.
        fprintf(output, "Points: ");
        if (!printCoords || window -> X == NULL) {
            fprintf(output, "null\n");
        } else {
            fprintf(output, "\n");
            for (int i = 0; i < N; i++) {
                print_row(output, window -> X[i], K, false);
                fprintf(output, "\n");
            }
        }
    }
}

// Make sure user defined arguments are valid.
// Accepts:
//  lodestarConfiguration_t* lodestar_config -> The configured parameters.
// Returns: int, 0 if all parameters are valid. -1 if user supplied an invalid value.
int check_configuration(LodestarConfiguration_t* lodestar_config) {
    // Check maf.
    if (lodestar_config -> maf < 0 || lodestar_config -> maf > 1) { 
        fprintf(stderr, "--maf %lf must be in [0, 1].\n", lodestar_config -> maf); 
        return -1;
    }
    // Check afMissing.
    if (lodestar_config -> afMissing < 0 || lodestar_config -> afMissing > 1) { 
        fprintf(stderr, "--afMissing %lf must be in [0, 1].\n", lodestar_config -> afMissing); 
        return -1;
    }
    // Check that input VCF file exists and region is valid if supplied.
    if (lodestar_config -> inputFileName == NULL) { 
        fprintf(stderr, "Must supply an input VCF file with -i.\n"); 
        return -1;
    }
    if (access(lodestar_config -> inputFileName, F_OK) != 0) {
        fprintf(stderr, "-i %s does not exist.\n", lodestar_config -> inputFileName);
        return -1;
    }
    // Check that output basename was given.
    if (lodestar_config -> outputBasename == NULL) { 
        fprintf(stderr, "Must supply an output basename with -o.\n"); 
        return -1;
    }
    // Check haplotype size is valid.
    if (lodestar_config -> HAP_SIZE <= 0) { 
        fprintf(stderr, "-h %d must be INT > 0.\n", lodestar_config -> HAP_SIZE); 
        return -1;
    }
    // Check window size is valid.
    if (!lodestar_config -> global && lodestar_config -> WINDOW_SIZE <= 0) { 
        fprintf(stderr, "-w %d must be INT > 0 supplied by the user.\n", lodestar_config -> WINDOW_SIZE); 
        return -1;
    }
    // Check step size is valid.
    if (!lodestar_config -> global && (lodestar_config -> STEP_SIZE <= 0 || lodestar_config -> STEP_SIZE > lodestar_config -> WINDOW_SIZE)) { 
        fprintf(stderr, "-s %d must be 0 < INT <= WINDOW_SIZE supplied by the user.\n", lodestar_config -> STEP_SIZE); 
        return -1;
    }
    // Check number of threads.
    if (lodestar_config -> threads <= 0) { 
        fprintf(stderr, "--threads %d must be INT > 0.\n", lodestar_config -> threads); 
        return -1;
    }
    // Check MAX_GAP.
    if (lodestar_config -> MAX_GAP < 1) { 
        fprintf(stderr, "--gap %d must be INT > 0.\n", lodestar_config -> MAX_GAP); 
        return -1;
    }
    // Check pthresh.
    if (lodestar_config -> pthresh < 0 || lodestar_config -> pthresh > 1) { 
        fprintf(stderr, "--pthresh %lf must be in (0, 1].\n", lodestar_config -> pthresh); 
        return -1;
    }
    // Check number of permutation.
    if (lodestar_config -> NUM_PERMS <= 0) { 
        fprintf(stderr, "--perms %d must be INT > 0.\n", lodestar_config -> NUM_PERMS); 
        return -1;
    }
    // Check tthresh.
    if (lodestar_config -> tthresh < 0 || lodestar_config -> tthresh > 1) { 
        fprintf(stderr, "--tthresh %lf must be in [0, 1].\n", lodestar_config -> tthresh); 
        return -1;
    }
    // Check target file if exists.
    if (lodestar_config -> targetFileName != NULL && access(lodestar_config -> targetFileName, F_OK) != 0) {
        fprintf(stderr, "--target %s does not exist.\n", lodestar_config -> targetFileName); 
    }
    // All parameters are valid, return success.
    return 0;
}

// Print help menu for LODESTAR.
// Accepts: void.
// Returns: void.
void print_help() {
    fprintf(stderr, "\n");
    fprintf(stderr, "LODESTAR v1.0\n");
    fprintf(stderr, "----------------------\n\n");
    fprintf(stderr, "Written by T. Quinn Smith\n");
    fprintf(stderr, "Principal Investigator: Zachary A. Szpiech\n");
    fprintf(stderr, "The Pennsylvania State University\n\n");
    fprintf(stderr, "Usage: lodestar [options] -w <WINDOW_SIZE> -s <STEP_SIZE> -i <input.vcf.gz> -o <outputBasename>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   --help                  Prints help menu and exits.\n");
    fprintf(stderr, "   --version               Prints version number and exits.\n");
    fprintf(stderr, "   -i STR                  Path to input vcf file.\n");
    fprintf(stderr, "   -o STR                  Output base names for files.\n");
    fprintf(stderr, "   -h INT                  Number of loci in a haplotype.\n");
    fprintf(stderr, "                               Default 1. Not used when --global is set.\n");
    fprintf(stderr, "   -w INT                  Number of haplotypes in a window.\n");
    fprintf(stderr, "                               Must be set by user. Not used when --global is set.\n");
    fprintf(stderr, "   -s INT                  Number of haplotypes to increment the sliding window.\n");
    fprintf(stderr, "                               Must be set by user. Not used when --global is set.\n");
    fprintf(stderr, "   -k INT                  Dimension to project samples into. Must be less than number of samples.\n");
    fprintf(stderr, "                               Default 2. Must be less than number of samples in VCF.\n");
    fprintf(stderr, "   --threads INT           Number of threads to use in computation.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   --similarity            Compute similarity between sets of points instead of dissimilarity.\n");
    fprintf(stderr, "   --global                Compute only the global set of points.\n");
    fprintf(stderr, "                               Takes presedence over windowing parameters.\n");
    fprintf(stderr, "   --target STR            A n-by-k csv file containing user defined coordinates to perform Procrustes analysis.\n");
    fprintf(stderr, "   --pthresh DOUBLE        Print coordinates of all windows less than or equal to threshold.\n");
    fprintf(stderr, "                               Window must also satisfy --tthresh. Default 0.\n");
    fprintf(stderr, "   --perms INT             The number of permutations to execute.\n");
    fprintf(stderr, "                               Default 10000. Permutation test does not execute if --pthresh is 0.\n");
    fprintf(stderr, "   --tthresh DOUBLE        Print coordinates of all windows greater than or equal to threshold.\n");
    fprintf(stderr, "                               Window must also satisfy --pthresh. Default 0.95.\n");
    fprintf(stderr, "   --maf DOUBLE            Drops biallelic VCF records with a MAF less than threshold.\n");
    fprintf(stderr, "                               Default 0.05\n");
    fprintf(stderr, "   --afMissing DOUBLE      Drops VCF records with fraction of genotypes missing greater than threshold.\n");
    fprintf(stderr, "                               Default 0.1\n");
    fprintf(stderr, "   --gap INT               Drops window if window covers more than gap value in basepairs.\n");
    fprintf(stderr, "                               Default 1000000\n");
    fprintf(stderr, "   --json                  Prints window information in JSON format instead of TXT.\n");
    fprintf(stderr, "\n");
}

// Long options used for LODESTAR.
static ko_longopt_t long_options[] = {
    {"help",            ko_no_argument,         300},
    {"version",         ko_no_argument,         'v'},
    {"similarity",      ko_no_argument,         302},
    {"global",          ko_no_argument,         303},
    {"threads",         ko_required_argument,   304},
    {"pthresh",         ko_required_argument,   305},
    {"perms",           ko_required_argument,   306},
    {"tthresh",         ko_required_argument,   307},
    {"maf",             ko_required_argument,   309},
    {"afMissing",       ko_required_argument,   310},
    {"json",            ko_no_argument,         314},
    {"target",          ko_required_argument,   315},
    {"gap",             ko_required_argument,   316},
    {"input",           ko_required_argument,   'i'},
    {"output",          ko_required_argument,   'o'},
    {"dimension",       ko_required_argument,   'k'},
    {"haplotype",       ko_required_argument,   'h'},
    {"window",          ko_required_argument,   'w'},
    {"step",            ko_required_argument,   's'},
    {0, 0, 0}
};

LodestarConfiguration_t* init_lodestar_config(int argc, char *argv[]) {

    // Print help menu when no options/arguments were given.
    if (argc == 1) {
        print_help();
        return NULL;
    } 

    // Parse user options.
    const char *opt_str = "i:o:h:w:s:k:v";
    ketopt_t options = KETOPT_INIT;
    int c;

    // Pass through options. Check for options requiring an argument that were not given one.
    //  Also, check if any options are supplied that are not defined. If user supplies the help
    //  argument, print help menu and exit, likewise for version.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return NULL;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return NULL;
            case 300: print_help(); return NULL;
            case 'v': fprintf(stderr, "Version 1.0 May 2024.\n"); return NULL;
        }
	}

    // Set defaults for LODESTAR configuration.
    LodestarConfiguration_t* lodestar_config = calloc(1, sizeof(LodestarConfiguration_t));
    lodestar_config -> inputFileName = NULL;
    lodestar_config -> outputBasename = NULL;
    lodestar_config -> HAP_SIZE = 1;
    lodestar_config -> WINDOW_SIZE = -1;
    lodestar_config -> STEP_SIZE = -1;
    lodestar_config -> k = 2;
    lodestar_config -> threads = 1;
    lodestar_config -> similarity = false;
    lodestar_config -> global = false;
    lodestar_config -> targetFileName = NULL;
    lodestar_config -> pthresh = 0;
    lodestar_config -> NUM_PERMS = 10000;
    lodestar_config -> tthresh = 0.95;
    lodestar_config -> maf = 0.05;
    lodestar_config -> afMissing = 0.1;
    lodestar_config -> MAX_GAP = 1000000;
    lodestar_config -> useJsonOutput = false;

    // Parse command line arguments.
    options = KETOPT_INIT;
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'i': lodestar_config -> inputFileName = options.arg; break;
            case 'o': lodestar_config -> outputBasename = options.arg; break;
            case 'h': lodestar_config -> HAP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'w': lodestar_config -> WINDOW_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 's': lodestar_config -> STEP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'k': lodestar_config -> k = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 302: lodestar_config -> similarity = true; break;
            case 303: lodestar_config -> global = true; break;
            case 304: lodestar_config -> threads = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 305: lodestar_config -> pthresh = strtod(options.arg, (char**) NULL); break;
            case 306: lodestar_config -> NUM_PERMS = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 307: lodestar_config -> tthresh = strtod(options.arg, (char**) NULL); break;
            case 309: lodestar_config -> maf = strtod(options.arg, (char**) NULL); break;
            case 310: lodestar_config -> afMissing = strtod(options.arg, (char**) NULL); break;
            case 314: lodestar_config -> useJsonOutput = true; break;
            case 315: lodestar_config -> targetFileName = options.arg; break;
            case 316: lodestar_config -> MAX_GAP = (int) strtol(options.arg, (char**) NULL, 10); break;
        }
	}

    // Check to that user input is all valid.
    if (check_configuration(lodestar_config) != 0)
        return NULL;

    return lodestar_config;

}