
// File: main.c
// Date: 9 May 2024
// Version 1: 17 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Command-line interface for LODESTAR and main analysis.

#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "SlidingWindow.h"
#include "ProcrustesAnalysis.h"
#include "Logger.h"
#include "../lib/ketopt.h"
#include "Matrix.h"
MATRIX_INIT(double, double)

// Used to index the lower triangle of a packed matrix.
#define INDEX(i, j, N) (i <= j ? i + j * (j + 1) / 2 : j + i * (i + 1) / 2)

// Print the single entry fields of a window.
// Accepts:
//  FILE* output -> The output stream to print to.
//  Window_t* window -> The window whose fields will be printed.
// Returns: void.
void print_window_info(FILE* output, Window_t* window) {
    fprintf(output, "%d\t", window -> winNum);
    fprintf(output, "%d\t", window -> winNumOnChrom);
    fprintf(output, "%s\t", ks_str(window -> chromosome));
    fprintf(output, "%d\t", window -> startCoord);
    fprintf(output, "%d\t", window -> endCoord);
    fprintf(output, "%d\t", window -> numLoci);
    fprintf(output, "%d\t", window -> numHaps);
    fprintf(output, "%lf\t", window -> pval);
    fprintf(output, "%lf\n", window -> t);
}

// Print an array.
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

// Print the attributes of a window.
//  Do not like how this method is written, but I do not think there is a way around it.
// Accepts:
//  FILE* output -> The output stream to print to.
//  kstring_t** sampleNames -> Array of sample names corresponding to each point in the window.
//              Only used when useLongOutput is set.
//  Window_t* window -> The window whose attributes will be printed.
//  int N -> The number of samples.
//  int K -> The dimension of the projected points.
//  bool useJsonOutput -> Print the attributes in JSON format. Otherwise, print in plain text.
//  bool useLongOutput -> Print coordinates and pairwise computations in long format.
//  bool asdToIbs -> Print pairwise computations as IBS counts. Otherwise, ASD.
//  bool printCoords -> Print the coordinates of the window if they exist.
// Returns: void.
void print_window_coords(FILE* output, kstring_t** sampleNames, Window_t* window, int N, int K, bool useJsonOutput, bool useLongOutput, bool asdToIbs, bool printCoords) {
    // Print JSON.
    if (useJsonOutput) {
        // Print single fields.
        fprintf(output, "\t{\n");
        fprintf(output, "\t\t\"Window Number\": %d,\n", window -> winNum);
        fprintf(output, "\t\t\"Window Number on Chromosome\": %d,\n", window -> winNumOnChrom);
        fprintf(output, "\t\t\"Chromosome\": %s,\n", ks_str(window -> chromosome));
        fprintf(output, "\t\t\"Start Coordinate\": %d,\n", window -> startCoord);
        fprintf(output, "\t\t\"End Coordinate\": %d,\n", window -> endCoord);
        fprintf(output, "\t\t\"Number of Loci\": %d,\n", window -> numLoci);
        fprintf(output, "\t\t\"Number of Haplotypes\": %d,\n", window -> numHaps);
        if (window -> pval == 0)
            fprintf(output, "\t\t\"P-Value\": null,\n");
        else 
            fprintf(output, "\t\t\"P-Value\": %lf,\n", window -> pval);
        if (window -> t == -1)
            fprintf(output, "\t\t\"t-statistic\": null,\n");
        else 
            fprintf(output, "\t\t\"t-statistic\": %lf,\n", window -> t);
        // Print points.
        fprintf(output, "\t\t\"Points\": ");
        if (!printCoords || window -> X == NULL) {
            fprintf(output, "null\n");
        } else {
            fprintf(output, "[\n");
            for (int i = 0; i < N; i++) {
                fprintf(output, "\t\t\t[");
                if (useLongOutput)
                    fprintf(output, "\"%s\", ", ks_str(sampleNames[i]));
                print_row(output, window -> X[i], K, true);
                if (i != N - 1)
                    fprintf(output, ",\n");
                else
                    fprintf(output, "\n");
            }
            fprintf(output, "\t\t]\n");
        }
        // Print pairwise.
        fprintf(output, "\t\t\"Pairwise\": ");
        if (!window -> saveIBS) {
            fprintf(output, "null\n");
        } else {
            fprintf(output, "[\n");
            if (useLongOutput) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < i; j++) {
                        fprintf(output, "\t\t\t[\"%s\", ", ks_str(sampleNames[i]));
                        fprintf(output, "\"%s\", ", ks_str(sampleNames[j]));
                        if (asdToIbs) 
                            fprintf(output, "%lf]", ibs_to_asd(window -> ibs[INDEX(i, j, N)]));
                        else
                            fprintf(output, "%d/%d/%d]", window -> ibs[INDEX(i, j, N)].ibs0, window -> ibs[INDEX(i, j, N)].ibs1, window -> ibs[INDEX(i, j, N)].ibs2);
                        if (i * j != (N - 1) * (N - 1))
                            fprintf(output, ",");
                        fprintf(output, "\n");
                    }
                }
                fprintf(output, "\t\t]\n");
            } else {
                fprintf(output, "[\n");
                for (int i = 0; i < N; i++) {
                    fprintf(output, "\t\t\t[");
                    for (int j = 0; j < i; j++) {
                        if (!asdToIbs)
                            fprintf(output, "%lf", ibs_to_asd(window -> ibs[INDEX(i, j, N)]));
                        else
                            fprintf(output, "%d/%d/%d", window -> ibs[INDEX(i, j, N)].ibs0, window -> ibs[INDEX(i, j, N)].ibs1, window -> ibs[INDEX(i, j, N)].ibs2);
                        if (j != i)
                            fprintf(output, ", ");
                    }
                    fprintf(output, "]");
                    if (i != N - 1)
                        fprintf(output, ",\n");
                    else 
                        fprintf(output, "\n");
                }
                fprintf(output, "\t\t]\n");
            }
        }
        fprintf(output, "\t}");
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
        if (window -> pval == 0)
            fprintf(output, "P-Value: null\n");
        else 
            fprintf(output, "P-Value: %lf\n", window -> pval);
        if (window -> t == -1)
            fprintf(output, "t-statistic: null\n");
        else
            fprintf(output, "t-statistic: %lf\n", window -> t);
        // Print points.
        fprintf(output, "Points: ");
        if (!printCoords || window -> X == NULL) {
            fprintf(output, "null\n");
        } else {
            fprintf(output, "\n");
            for (int i = 0; i < N; i++) {
                if (useLongOutput)
                    fprintf(output, "%s\t", ks_str(sampleNames[i]));
                print_row(output, window -> X[i], K, false);
                fprintf(output, "\n");
            }
        }
        // Print pairwise.
        fprintf(output, "Pairwise: ");
        if (!window -> saveIBS) {
            fprintf(output, "null\n");
        } else {
            if (useLongOutput) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < i; j++) {
                        fprintf(output, "%s\t", ks_str(sampleNames[i]));
                        fprintf(output, "%s\t", ks_str(sampleNames[j]));
                        if (asdToIbs)
                            fprintf(output, "%lf\n", ibs_to_asd(window -> ibs[INDEX(i, j, N)]));
                        else
                            fprintf(output, "%d/%d/%d\n", window -> ibs[INDEX(i, j, N)].ibs0, window -> ibs[INDEX(i, j, N)].ibs1, window -> ibs[INDEX(i, j, N)].ibs2);
                    }
                }
            } else {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < i; j++) {
                        if (!asdToIbs)
                            fprintf(output, "%lf", ibs_to_asd(window -> ibs[INDEX(i, j, N)]));
                        else
                            fprintf(output, "%d/%d/%d", window -> ibs[INDEX(i, j, N)].ibs0, window -> ibs[INDEX(i, j, N)].ibs1, window -> ibs[INDEX(i, j, N)].ibs2);
                        if (j != i)
                            fprintf(output, "\t");
                    }
                    fprintf(output, "\n");
                }
            }
        }
    }
}

// Defines all the command line options supplied
//  by user and parameters to run the LODESTAR analysis.
typedef struct {
    // The input VCF file.
    char* inputFileName;
    // The parser object used to read the VCF file.
    VCFLocusParser_t* parser;
    // The basename for output files.
    char* outputBasename;
    // The haplotype size in number of loci.
    int HAP_SIZE;
    // The windows size in number of haplotypes.
    int WINDOW_SIZE;
    // The step size in number of haplotypes.
    int STEP_SIZE;
    // The dimension of the reduced points.
    int k;
    // The number of threads to use.
    int threads;
    // Flag to indicate if Procrustes analysis is measuring
    //  similarity or dissimilarity.
    bool similarity;
    // Flag to indicate if we are just calculating global window.
    bool global;
    // Name of file of user specified points to perform Procrustes
    //  analysis with.
    char* targetFileName;
    // Pointer to file of user specified points to perform Procrustes
    //  analysis with.
    FILE* targetFile;
    // Pvalue threshold used to print points.
    double pthresh;
    // Number of permutations to execute in permutation test.
    int NUM_PERMS;
    // t-statistic threshold used to print points.
    double tthresh;
    // Regions to include/exclude while parsing VCF file.
    char* regions;
    // Minor allele frequency threshhold.
    double maf;
    // Missing allele frequency threshold.
    double afMissing;
    // Regions to save window IBS/ASD matrix.
    char* saveRegionsStr;
    RegionSet_t* saveRegions;
    // Covnert ASD values to IBS counts.
    bool asdToIbs;
    // Regions to print coodinates of samples in dimension-k.
    char* printRegionsStr;
    RegionSet_t* printRegions;
    // Output long format instead of lower triangle for pairwise calculations.
    bool useLongOutput;
    // Output file in JSON format instead of text file.
    bool useJsonOutput;
} LodestarConfiguration_t;

// Make sure user defined arguments are valid.
// Accepts:
//  LodestarConfiguration_t* lodestarConfig -> The configured parameters.
// Returns: int, 0 if all parameters are valid. 1 if user supplied an invalid value.
int check_configuration(LodestarConfiguration_t* lodestarConfig) {
    // Check maf.
    if (lodestarConfig -> maf < 0 || lodestarConfig -> maf > 1) { fprintf(stderr, "--maf must be in [0, 1].\n"); return 1;}
    // Check afMissing.
    if (lodestarConfig -> afMissing < 0 || lodestarConfig -> afMissing > 1) { fprintf(stderr, "--afMissing must be in [0, 1].\n"); return 1;}
    // Check that input VCF file exists and region is valid if supplied.
    if (lodestarConfig -> inputFileName == NULL) { fprintf(stderr, "Must supply an input VCF file with -i.\n"); return 1;}
    kstring_t* r = init_kstring(NULL);
    bool takeComplement; 
    if (lodestarConfig -> regions != NULL) {
        takeComplement = lodestarConfig -> regions[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestarConfig -> regions + 1, r);
        else 
            ks_overwrite(lodestarConfig -> regions, r);
        if ((lodestarConfig -> parser = init_vcf_locus_parser(lodestarConfig -> inputFileName, r, takeComplement, lodestarConfig -> maf, lodestarConfig -> afMissing, true)) == NULL) {
            destroy_kstring(r);
            fprintf(stderr, "-i file does not exist or --regions invalid.\n");
            return 1;
        }
    } else {
        if ((lodestarConfig -> parser = init_vcf_locus_parser(lodestarConfig -> inputFileName, NULL, false, lodestarConfig -> maf, lodestarConfig -> afMissing, true)) == NULL) {
            destroy_kstring(r);
            fprintf(stderr, "-i file does not exist.\n");
            return 1;
        }
    }
    // Check regions to save is valid.
    if (lodestarConfig -> saveRegionsStr != NULL) {
        takeComplement = lodestarConfig -> saveRegionsStr[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestarConfig -> saveRegionsStr + 1, r);
        else 
            ks_overwrite(lodestarConfig -> saveRegionsStr, r);
        if ((lodestarConfig -> saveRegions = init_region_set(r, takeComplement)) == NULL) {
            destroy_kstring(r);
            fprintf(stderr, "--save region is invalid.\n");
            return 1;
        }
    }
    // Check region to print is valid.
    if (lodestarConfig -> printRegionsStr != NULL) {
        takeComplement = lodestarConfig -> printRegionsStr[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestarConfig -> printRegionsStr + 1, r);
        else 
            ks_overwrite(lodestarConfig -> printRegionsStr, r);
        if ((lodestarConfig -> printRegions = init_region_set(r, takeComplement)) == NULL) {
            destroy_kstring(r);
            fprintf(stderr, "--printCoords region is invalid.\n");
            return 1;
        }
    }
    destroy_kstring(r);
    // Check that output basename was given.
    if (lodestarConfig -> outputBasename == NULL) { fprintf(stderr, "Must supply an output basename with -o.\n"); return 1;}
    // Check haplotype size is valid.
    if (lodestarConfig -> HAP_SIZE <= 0) { fprintf(stderr, "-h must be INT > 0.\n"); return 1;}
    // Check window size is valid.
    if (!lodestarConfig -> global && lodestarConfig -> WINDOW_SIZE <= 0) { fprintf(stderr, "-w must be INT > 0 supplied by the user.\n"); return 1;}
    // Check step size is valid.
    if (!lodestarConfig -> global && (lodestarConfig -> STEP_SIZE <= 0 || lodestarConfig -> STEP_SIZE > lodestarConfig -> WINDOW_SIZE)) { fprintf(stderr, "-s must be 0 < INT <= WINDOW_SIZE supplied by the user.\n"); return 1;}
    // Check k.
    if (lodestarConfig -> k <= 0 || lodestarConfig -> k >= lodestarConfig -> parser -> numSamples) { fprintf(stderr, "-k must be 1 <= INT < N.\n"); return 1;}
    // Check number of threads.
    if (lodestarConfig -> threads <= 0) { fprintf(stderr, "--threads must be INT > 0.\n"); return 1;}
    // Check pthresh.
    if (lodestarConfig -> pthresh < 0 || lodestarConfig -> pthresh > 1) { fprintf(stderr, "--pthresh must be in (0, 1].\n"); return 1;}
    // Check number of permutation.
    if (lodestarConfig -> NUM_PERMS <= 0) { fprintf(stderr, "--perms must be INT > 0.\n"); return 1;}
    // Check tthresh.
    if (lodestarConfig -> tthresh < 0 || lodestarConfig -> tthresh > 1) { fprintf(stderr, "--tthresh must be in [0, 1].\n"); return 1;}
    // Check target file if exists.
    if (lodestarConfig -> targetFileName != NULL) {
        lodestarConfig -> targetFile = fopen(lodestarConfig -> targetFileName, "r");
        if (lodestarConfig -> targetFile == NULL) {
            fprintf(stderr, "--target %s does not exist.\n", lodestarConfig -> targetFileName);
            return 1;
        }
        int numLines = 0, dim = 1;
        char c;
        while ((c = getc(lodestarConfig -> targetFile)) != EOF) {
            if (c == '\t')
                dim++;
            if (c == '\n') {
                // Check dimension of target file.
                if (dim != lodestarConfig -> k) {
                    fprintf(stderr, "Dimension of --target %s does not match that of k.\n", lodestarConfig -> targetFileName);
                    fclose(lodestarConfig -> targetFile);
                    return 1;
                }
                dim = 1;
                numLines++;
            }
        }
        // Check the number of samples in target file.
        if (numLines != lodestarConfig -> parser -> numSamples) {
            fprintf(stderr, "Number of samples in --target %s does not match that of the input file.\n", lodestarConfig -> targetFileName);
            return 1;
        }
    }
    // All parameters are valid, return success.
    return 0;
}

#define PRINT_BOOL(X) (X ? "true" : "false")

// Print verbose LODESTAR configuration.
// Accepts:
//  FILE* output -> The output stream to print the configuration.
//  LodestarConfiguration_t lodestart_config -> THe configuration to print.
// Returns: void.
void print_configuration(FILE* output, LodestarConfiguration_t lodestarConfig) {
    fprintf(output, "LODESTAR Configuration\n");
    fprintf(output, "----------------------\n");
    fprintf(output, "Input file: %s\n", lodestarConfig.inputFileName);
    fprintf(output, "Output basename: %s\n", lodestarConfig.outputBasename);
    fprintf(output, "Haplotype size: %d\n", lodestarConfig.HAP_SIZE);
    if (!lodestarConfig.global) {
        fprintf(output, "Window size: %d\n", lodestarConfig.WINDOW_SIZE);
        fprintf(output, "Step size: %d\n", lodestarConfig.STEP_SIZE);
    }
    fprintf(output, "K: %d\n", lodestarConfig.k);
    fprintf(output, "Use similarity: %s\n", PRINT_BOOL(lodestarConfig.similarity));
    fprintf(output, "Calculate genome-wide only: %s\n", PRINT_BOOL(lodestarConfig.global));
    fprintf(output, "Number of threads used: %d\n", lodestarConfig.threads);
    fprintf(output, "P-value threshold (Disabled if 0): %lf\n", lodestarConfig.pthresh);
    fprintf(output, "Number of permutations (Disabled when p-value is 0): %d\n", lodestarConfig.NUM_PERMS);
    fprintf(output, "Procrustes statistic threshold: %lf\n", lodestarConfig.tthresh);
    if (lodestarConfig.regions != NULL)
        fprintf(output, "Include/exclude records from input: %s\n", lodestarConfig.regions);
    fprintf(output, "Biallelic minor allele frequency threshold: %lf\n", lodestarConfig.maf);
    fprintf(output, "Missing allele frequency threshold: %lf\n", lodestarConfig.afMissing);
    fprintf(output, "Convert ASD values to IBS counts: %s\n", PRINT_BOOL(lodestarConfig.asdToIbs));
    if (lodestarConfig.saveRegionsStr != NULL)
        fprintf(output, "Save IBS/ASD values for windows overlapping/not overlapping: %s\n", lodestarConfig.saveRegionsStr);
    if (lodestarConfig.printRegionsStr != NULL)
        fprintf(output, "Print coordinates of windows that are/are not overlapping: %s\n", lodestarConfig.printRegionsStr);
    if (lodestarConfig.targetFileName != NULL)
        fprintf(output, "File of coordinates to perform Procrustes analysis against: %s\n", lodestarConfig.targetFileName);
    fprintf(output, "Use long-format output: %s\n", PRINT_BOOL(lodestarConfig.useLongOutput));
    fprintf(output, "Save as JSON file: %s\n", PRINT_BOOL(lodestarConfig.useJsonOutput));
}

// Destroy all dynamically allocated memory associated with lodestarConfig.
// Accepts:
//  LodestarConfiguration_t lodestarConfig.The configuration to destroy.
// Returns: void.
void destroy_lodestarConfiguration(LodestarConfiguration_t lodestarConfig) {
    if (lodestarConfig.parser != NULL)
        destroy_vcf_locus_parser(lodestarConfig.parser);
    if (lodestarConfig.printRegions != NULL)
        destroy_region_set(lodestarConfig.printRegions);
    if (lodestarConfig.saveRegions != NULL)
        destroy_region_set(lodestarConfig.saveRegions);
}

// Print help menu for LODESTAR.
// Accepts: void.
// Returns: void.
void print_help() {
    fprintf(stderr, "\n");
    fprintf(stderr, "LODESTAR v1.0 May 2024\n");
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
    fprintf(stderr, "   --regions [^]REGS       Include/exclude records from VCF defined by REGS.\n");
    fprintf(stderr, "                               Default NULL.\n");
    fprintf(stderr, "   --printCoords [^]REGS   Print coordinates that are overlapping/non-overlapping with REGS.\n");
    fprintf(stderr, "                               Default NULL.\n");
    fprintf(stderr, "   --save [^]REGS          Save IBS/ASD values for windows that are overlapping/non-overlapping with REGS.\n");
    fprintf(stderr, "                               Default NULL. See --asdToIbs\n");
    fprintf(stderr, "   --asdToIbs              Convert IBS values to ASD values in output.\n");
    fprintf(stderr, "                               Default false.\n");
    fprintf(stderr, "   --maf DOUBLE            Drops biallelic VCF records with a MAF less than threshold.\n");
    fprintf(stderr, "                               Default 0.05\n");
    fprintf(stderr, "   --afMissing DOUBLE      Drops VCF records with fraction of genotypes missing greater than threshold.\n");
    fprintf(stderr, "                               Default 0.5.\n");
    fprintf(stderr, "   --long                  Prints calculations in long format instead of matrix form.\n");
    fprintf(stderr, "   --json                  Prints window information in JSON format instead of TXT.\n");
    fprintf(stderr, "Types:\n");
    fprintf(stderr, "   STR                     A string.\n");
    fprintf(stderr, "   INT                     A non-negative integer.\n");
    fprintf(stderr, "   DOUBLE                  A real number between 0 and 1, inclusive.\n");
    fprintf(stderr, "   REGS                    REG,REGS | REG\n");
    fprintf(stderr, "   REG                     STR | STR:INT | STR:-INT | STR:INT- | STR:INT-INT\n");
    fprintf(stderr, "   Note: ^ before REGS denotes taking the complement of the intervals over the whole genome.\n");
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
    {"regions",         ko_required_argument,   308},
    {"maf",             ko_required_argument,   309},
    {"afMissing",       ko_required_argument,   310},
    {"save",            ko_required_argument,   311},
    {"asdToIbs",        ko_no_argument,         312},
    {"long",            ko_no_argument,         313},
    {"json",            ko_no_argument,         314},
    {"target",          ko_required_argument,   315},
    {"printCoords",     ko_required_argument,   316},
    {"input",           ko_required_argument,   'i'},
    {"output",          ko_required_argument,   'o'},
    {"dimension",       ko_required_argument,   'k'},
    {"haplotype",       ko_required_argument,   'h'},
    {"window",          ko_required_argument,   'w'},
    {"step",            ko_required_argument,   's'},
    {0, 0, 0}
};

int main (int argc, char *argv[]) {

    // Seed random number generator.
    srand(time(NULL));

    // Single character aliases for long options.
    const char *opt_str = "i:o:h:w:s:k:v";
    ketopt_t options = KETOPT_INIT;
    int c;

    // Print help menu when no options/arguments were given.
    if (argc == 1) {
        print_help();
        return 1;
    }        

    // Pass through options. Check for options requiring an argument that were not given one.
    //  Also, check if any options are supplied that are not defined. If user supplies the help
    //  argument, print help menu and exit, likewise for version.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return 1;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return 1;
            case 300: print_help(); return 0;
            case 'v': fprintf(stderr, "Version 1.0 May 2024.\n"); return 0;
        }
	}
	options = KETOPT_INIT;

    // Set defaults for LODESTAR configuration.
    LodestarConfiguration_t lodestarConfig;
    lodestarConfig.inputFileName = NULL;
    lodestarConfig.parser = NULL;
    lodestarConfig.outputBasename = NULL;
    lodestarConfig.HAP_SIZE = 1;
    lodestarConfig.WINDOW_SIZE = -1;
    lodestarConfig.STEP_SIZE = -1;
    lodestarConfig.k = 2;
    lodestarConfig.threads = 1;
    lodestarConfig.similarity = false;
    lodestarConfig.global = false;
    lodestarConfig.targetFileName = NULL;
    lodestarConfig.targetFile = NULL;
    lodestarConfig.pthresh = 0;
    lodestarConfig.NUM_PERMS = 10000;
    lodestarConfig.tthresh = 0.95;
    lodestarConfig.regions = NULL;
    lodestarConfig.maf = 0.05;
    lodestarConfig.afMissing = 0.5;
    lodestarConfig.saveRegionsStr = NULL;
    lodestarConfig.saveRegions = NULL;
    lodestarConfig.asdToIbs = false;
    lodestarConfig.printRegionsStr = NULL;
    lodestarConfig.printRegions = NULL;
    lodestarConfig.useLongOutput = false;
    lodestarConfig.useJsonOutput = false;

    // Parse command line arguments.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'i': lodestarConfig.inputFileName = options.arg; break;
            case 'o': lodestarConfig.outputBasename = options.arg; break;
            case 'h': lodestarConfig.HAP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'w': lodestarConfig.WINDOW_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 's': lodestarConfig.STEP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'k': lodestarConfig.k = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 302: lodestarConfig.similarity = true; break;
            case 303: lodestarConfig.global = true; break;
            case 304: lodestarConfig.threads = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 305: lodestarConfig.pthresh = strtod(options.arg, (char**) NULL); break;
            case 306: lodestarConfig.NUM_PERMS = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 307: lodestarConfig.tthresh = strtod(options.arg, (char**) NULL); break;
            case 308: lodestarConfig.regions = options.arg; break;
            case 309: lodestarConfig.maf = strtod(options.arg, (char**) NULL); break;
            case 310: lodestarConfig.afMissing = strtod(options.arg, (char**) NULL); break;
            case 311: lodestarConfig.saveRegionsStr = options.arg; break;
            case 312: lodestarConfig.asdToIbs = true; break;
            case 313: lodestarConfig.useLongOutput = true; break;
            case 314: lodestarConfig.useJsonOutput = true; break;
            case 315: lodestarConfig.targetFileName = options.arg; break;
            case 316: lodestarConfig.printRegionsStr = options.arg; break;
        }
	}
    
    // Check configuration. If invalid argument, exit program.
    if (check_configuration(&lodestarConfig) != 0) {
        fprintf(stderr, "Exiting!\n");
        return 1;
    }

    // If this point is reached, the configuration is valid, and
    //  we can execute our analysis.

    // Setup the logfile.
    kstring_t* outputBasename = init_kstring(lodestarConfig.outputBasename);
    kputs(".log", outputBasename);
    INIT_LOG(ks_str(outputBasename));
    printf("\nLogging progress in %s\n\n", ks_str(outputBasename));
    ks_overwrite(lodestarConfig.outputBasename, outputBasename);
    // Print configuration to log file.
    print_configuration(logger -> file, lodestarConfig);
    // Print sample names in log file.
    fprintf(logger -> file, "\nSample Names:\n");
    for (int i = 0; i < lodestarConfig.parser -> numSamples; i++)
        fprintf(logger -> file, "%s\n", ks_str(lodestarConfig.parser -> sampleNames[i]));
    fprintf(logger -> file, "\n");

    // Create the haplotype encoder.
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(lodestarConfig.parser -> numSamples);

    Window_t* global = NULL;
    Window_t** windows = NULL;
    int numWindows = 0;

    // If we are only calculating the genome-wide ...
    if (lodestarConfig.global) {
        LOG_INFO("Performing genome-wide calculations ...\n");
        printf("Performing genome-wide calculations ...\n\n");

        global = global_window(lodestarConfig.parser, encoder, lodestarConfig.k, lodestarConfig.HAP_SIZE, lodestarConfig.threads);

        LOG_INFO("Finished genome-wide calculations ...\n");
        printf("Finished genome-wide calculations ...\n\n");
    // If we are calculating the sliding window ...
    } else {
        LOG_INFO("Performing sliding-window calculations ...\n");
        printf("Performing sliding-window calculations ...\n\n");

        windows = sliding_window(lodestarConfig.parser, encoder, lodestarConfig.saveRegions, lodestarConfig.k, lodestarConfig.HAP_SIZE, lodestarConfig.STEP_SIZE, lodestarConfig.WINDOW_SIZE, lodestarConfig.threads, &numWindows);

        // Only print pairwise for global if we are printing pairwise for other windows.
        if (lodestarConfig.saveRegionsStr == NULL) {
            free(windows[0] -> ibs);
            windows[0] -> ibs = NULL;
            windows[0] -> saveIBS = false;
        }

        LOG_INFO("Finished sliding-window calculations ...\n");
        printf("Finished sliding-window calculations ...\n\n");
    }

    double** target = NULL;
    double* target0 = NULL;

    // If a target file is supplied, open the file and create the target set of points.
    if (lodestarConfig.targetFileName != NULL) {
        lodestarConfig.targetFile = fopen(lodestarConfig.targetFileName, "r");
        target = create_matrix(double, encoder -> numSamples, lodestarConfig.k);
        target0 = (double*) calloc(lodestarConfig.k, sizeof(double));
        double value;
        for (int i = 0; i < encoder -> numSamples; i++) {
            for (int j = 0; j < lodestarConfig.k; j++) {
                fscanf(lodestarConfig.targetFile, "%lf", &value);
                target[i][j] = value;
            }
        }
        // Center user defined points.
        center_matrix(target, target0, encoder -> numSamples, lodestarConfig.k);
    } 
    if (lodestarConfig.targetFileName == NULL && windows != NULL) {
        target = windows[0] -> X;
    }

    // Used to transform points.
    RealSymEigen_t* eigen = init_real_sym_eigen(encoder -> numSamples);

    if ((lodestarConfig.global && global -> X != NULL && target != NULL) || target != NULL) {
        LOG_INFO("Beginning Procrustes Analysis ...\n");
        printf("Beginning Procrustes Analysis ...\n\n");
        // If only global was calculated and the user suppplied coordinates, perform Procrustes against the two.
        if (lodestarConfig.global) {
            // If we are executing Procrustes analysis once, do not multithread even if user allowed multithreading.
            double t = procrustes_statistic(global -> X, NULL, target, NULL, eigen, eigen -> N, lodestarConfig.k, false, lodestarConfig.similarity);
            global -> t = t;
            // Execute permutation test if user set p-value threshold.
            if (lodestarConfig.pthresh != 0) {
                double** shuffleX = create_matrix(double, encoder -> numSamples, lodestarConfig.k);
                global -> pval = permutation_test(global -> X, target, shuffleX, eigen, eigen -> N, lodestarConfig.k, lodestarConfig.similarity, t, lodestarConfig.NUM_PERMS);
                destroy_matrix(double, shuffleX, encoder -> numSamples);
            }
        } 
        
        procrustes_sliding_window(windows, numWindows, target, encoder -> numSamples, lodestarConfig.k, lodestarConfig.similarity, lodestarConfig.pthresh == 0 ? 0 : lodestarConfig.NUM_PERMS);
        
        LOG_INFO("Finished Procrustes Analysis ...\n");
        printf("Finished Procrustes Analysis ...\n\n");
    }

    // Error messages to signal that Procrustes analysis was not executed.
    if (lodestarConfig.global && global -> X == NULL && target != NULL) {
        LOG_ERROR("Could not perform Procrustes Analysis between genome-wide and user defined coordinates because genome-wide coordinates were of low rank.\n");
    } else if (!lodestarConfig.global && target == NULL) {
        LOG_ERROR("Could not perform Procrustes Analysis between windows and genome-wide coordinates because genome-wide coordinates were of low rank.\n");
    }

    LOG_INFO("Saving results to output files ...\n");
    printf("Saving results to output files ...\n\n");

    // Our output files.
    FILE* windowSummaries = NULL;
    FILE* windowCoords = NULL;
    FILE* globalCoords = NULL;
    
    // If we performed a sliding window.
    if (!lodestarConfig.global) {
        // Open long window output file.
        kputs("_windows", outputBasename);
        kputs(lodestarConfig.outputBasename, outputBasename);
        if (lodestarConfig.useJsonOutput)
            kputs(".json", outputBasename);
        else
            kputs(".txt", outputBasename);
        windowCoords = fopen(ks_str(outputBasename), "w");
        // Open summary file and print headers to summary file.
        kputs("_summary", outputBasename);
        kputs(lodestarConfig.outputBasename, outputBasename);
        kputs(".tsv", outputBasename);
        windowSummaries = fopen(ks_str(outputBasename), "w");
        fprintf(windowSummaries, "#Command: ");
        for (int i = 0; i < argc; i++) 
            fprintf(windowSummaries, "%s ", argv[i]);
        fprintf(windowSummaries, "\n");
        fprintf(windowSummaries, "Win\t");
        fprintf(windowSummaries, "WinChr\t");
        fprintf(windowSummaries, "Chr\t");
        fprintf(windowSummaries, "Start\t");
        fprintf(windowSummaries, "End\t");
        fprintf(windowSummaries, "nLoci\t");
        fprintf(windowSummaries, "nHaps\t");
        fprintf(windowSummaries, "p-val\t");
        fprintf(windowSummaries, "t-stat\n");
        if (lodestarConfig.useJsonOutput)
            fprintf(windowCoords, "[\n");
        bool printCoords;
        // For each window ...
        for (int i = 1; i < numWindows; i++) {
            // Print summary.
            print_window_info(windowSummaries, windows[i]);
            if (lodestarConfig.useJsonOutput)
                fprintf(windowCoords, ",");
            fprintf(windowCoords, "\n");
            // Determine if coordinates for the file should be printed and print to windows file.
            printCoords = (windows[i] -> X != NULL) && ((lodestarConfig.pthresh != 0 && windows[i] -> pval < lodestarConfig.pthresh) || (windows[i] -> t >= lodestarConfig.tthresh) || (lodestarConfig.printRegions != NULL && query_overlap(lodestarConfig.printRegions, windows[i] -> chromosome, windows[i] -> startCoord, windows[i] -> endCoord)));
            // Transform coordinates if we are printing points.
            if (printCoords)
                procrustes_statistic(windows[i] -> X, NULL, target, target0, eigen, encoder -> numSamples, lodestarConfig.k, true, lodestarConfig.similarity);
            print_window_coords(windowCoords, lodestarConfig.parser -> sampleNames, windows[i], encoder -> numSamples, lodestarConfig.k, lodestarConfig.useJsonOutput, lodestarConfig.useLongOutput, lodestarConfig.asdToIbs, printCoords);
        }
        if (lodestarConfig.useJsonOutput)
            fprintf(windowCoords, "\n]\n");
    }

    // Create global output file.
    kputs("_global", outputBasename);
    kputs(lodestarConfig.outputBasename, outputBasename);
    if (lodestarConfig.useJsonOutput)
        kputs(".json", outputBasename);
    else
        kputs(".txt", outputBasename);
    globalCoords = fopen(ks_str(outputBasename), "w");

    if (lodestarConfig.global && global -> X != NULL) {
        if (target != NULL)
            procrustes_statistic(global -> X, NULL, target, target0, eigen, encoder -> numSamples, lodestarConfig.k, true, lodestarConfig.similarity);
        print_window_coords(globalCoords, lodestarConfig.parser -> sampleNames, global, encoder -> numSamples, lodestarConfig.k, lodestarConfig.useJsonOutput, lodestarConfig.useLongOutput, lodestarConfig.asdToIbs, true);
    } 
    if (!lodestarConfig.global && windows[0] -> X != NULL) {
        if (target != NULL && target != windows[0] -> X) 
            procrustes_statistic(windows[0] -> X, NULL, target, target0, eigen, encoder -> numSamples, lodestarConfig.k, true, lodestarConfig.similarity);
        print_window_coords(globalCoords, lodestarConfig.parser -> sampleNames, windows[0], encoder -> numSamples, lodestarConfig.k, lodestarConfig.useJsonOutput, lodestarConfig.useLongOutput, lodestarConfig.asdToIbs, true);
    }

    LOG_INFO("Finished Analysis! Exiting ...\n");
    printf("Finished Analysis! Exiting ...\n\n");

    // Close all files and free all memory used in analysis.
    CLOSE_LOG();
    destroy_real_sym_eigen(eigen);
    destroy_kstring(outputBasename);
    if (lodestarConfig.targetFileName != NULL) {
        destroy_matrix(double, target, encoder -> numSamples);
        free(target0);
    }
    destroy_lodestarConfiguration(lodestarConfig);
    if (windowSummaries != NULL)
        fclose(windowSummaries);
    if (windowCoords != NULL)
        fclose(windowCoords);
    if (globalCoords != NULL)
        fclose(globalCoords);
    if (global != NULL)
        destroy_window(global, encoder -> numSamples);
    if (windows != NULL) {
        for (int i = 0; i < numWindows; i++)
            destroy_window(windows[i], encoder -> numSamples);
        free(windows);
    }
    destroy_haplotype_encoder(encoder);
    return 0;

}