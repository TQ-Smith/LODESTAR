
// File: main.c
// Date: 9 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Command-line interface for LODESTAR and main analysis.

#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "SlidingWindow.h"
#include "ProcrustesAnalysis.h"
#include "Logger.h"
#include "../lib/ketopt.h"

#define INDEX(i, j, N) (i <= j ? i + j * (j + 1) / 2 : j + i * (i + 1) / 2)

void print_window_info(FILE* output, Window_t* window) {
    fprintf(output, "%d\t", window -> winNum);
    fprintf(output, "%d\t", window -> winNumOnChrom);
    fprintf(output, "%s\t", ks_str(window -> chromosome));
    fprintf(output, "%d\t", window -> startCoord);
    fprintf(output, "%d\t", window -> endCoord);
    fprintf(output, "%d\t", window -> numLoci);
    if (window -> pval == 0)
        fprintf(output, "NA\t");
    else 
        fprintf(output, "%lf\t", window -> pval);
    fprintf(output, "%lf\n", window -> t);
}

void print_row(FILE* output, double* row, int K, bool json_output) {
    if (json_output) {
        for (int i = 0; i < K - 1; i++) {
            fprintf(output, "%lf, ", row[i]);
        }
        fprintf(output, "%lf]", row[K - 1]);
    } else {
        for (int i = 0; i < K - 1; i++) {
            fprintf(output, "%lf\t", row[i]);
        }
        fprintf(output, "%lf", row[K - 1]);
    }
}

void print_window_coords(FILE* output, kstring_t* sampleNames, Window_t* window, int N, int K, bool json_output, bool long_output, bool asdToIbs) {
    if (json_output) {
        fprintf(output, "{\n");
        fprintf(output, "\t\"Window Number\": %d,\n", window -> winNum);
        fprintf(output, "\t\"Window Number on Chromosome\": %d,\n", window -> winNumOnChrom);
        fprintf(output, "\t\"Chromosome\": %s,\n", ks_str(window -> chromosome));
        fprintf(output, "\t\"Start Coordinate\": %d,\n", window -> startCoord);
        fprintf(output, "\t\"End Coordinate\": %d,\n", window -> endCoord);
        fprintf(output, "\t\"Number of Loci\": %d,\n", window -> numLoci);
        if (window -> pval == 0)
            fprintf(output, "\t\"P-Value\": \"NA\",\n");
        else 
            fprintf(output, "\t\"P-Value\": %lf,\n", window -> pval);
        fprintf(output, "\t\"t-statistic\": %lf,\n", window -> t);
        fprintf(output, "\t\"Points\": ");
        if (window -> X == NULL) {
            fprintf(output, "\"NA\"\n");
        } else {
            fprintf(output, "[\n");
            for (int i = 0; i < N; i++) {
                fprintf(output, "\t\t[");
                if (long_output)
                    fprintf(output, "%s, ", sampleNames[i].s);
                print_row(output, window -> X[i], K, true);
                if (i != N - 1)
                    fprintf(output, ",");
            }
        }
        fprintf(output, "\t\"Pairwise\": ");
        if (!window -> saveIBS) {
            fprintf(output, "\"NA\"\n");
        } else {
            if (long_output) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j <= i; j++) {
                        fprintf(output, "%s, ", sampleNames[i].s);
                        fprintf(output, "%s, ", sampleNames[j].s);
                        if (asdToIbs)
                            fprintf(output, "%lf]", ibs_to_asd(window -> ibs[INDEX(i, j, N)]));
                        else
                            fprintf(output, "%d/%d/%d]", window -> ibs[INDEX(i, j, N)].ibs0, window -> ibs[INDEX(i, j, N)].ibs1, window -> ibs[INDEX(i, j, N)].ibs2);
                    }
                    if (i != N - 1)
                        fprintf(output, ",");
                    fprintf(output, "\n");
                }
            } else {
                for (int i = 0; i < N; i++) {
                    fprintf(output, "[");
                    for (int j = 0; j <= i; j++) {
                        if (asdToIbs)
                            fprintf(output, "%lf", ibs_to_asd(window -> ibs[INDEX(i, j, N)]));
                        else
                            fprintf(output, "%d/%d/%d", window -> ibs[INDEX(i, j, N)].ibs0, window -> ibs[INDEX(i, j, N)].ibs1, window -> ibs[INDEX(i, j, N)].ibs2);
                        if (j != i)
                            fprintf(output, ",");
                    }
                    fprintf(output, "]");
                    if (i != N - 1)
                        fprintf(output, ",\n");
                    else 
                        fprintf(output, "\n");
                }
            }
        }
        fprintf(output, "}\n");
    } else {
        fprintf(output, "Window Number: %d\n", window -> winNum);
        fprintf(output, "Window Number on Chromosome: %d\n", window -> winNumOnChrom);
        fprintf(output, "Chromosome: %s\n", ks_str(window -> chromosome));
        fprintf(output, "Start Coordinate: %d\n", window -> startCoord);
        fprintf(output, "End Coordinate: %d\n", window -> endCoord);
        fprintf(output, "Number of Loci: %d\n", window -> numLoci);
        if (window -> pval == 0)
            fprintf(output, "P-Value: NA\n");
        else 
            fprintf(output, "P-Value: %lf\n", window -> pval);
        fprintf(output, "t-statistic: %lf\n", window -> t);
        if (window -> X == NULL) {
            fprintf(output, "NA\n");
        } else {
            for (int i = 0; i < N; i++) {
                if (long_output)
                    fprintf(output, "%s\t", sampleNames[i].s);
                print_row(output, window -> X[i], K, false);
                fprintf(output, "\n");
            }
        }
        fprintf(output, "Pairwise: ");
        if (!window -> saveIBS) {
            fprintf(output, "NA\n");
        } else {
            if (long_output) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j <= i; j++) {
                        fprintf(output, "%s\t", sampleNames[i].s);
                        fprintf(output, "%s\t", sampleNames[j].s);
                        if (asdToIbs)
                            fprintf(output, "%lf\n", ibs_to_asd(window -> ibs[INDEX(i, j, N)]));
                        else
                            fprintf(output, "%d/%d/%d\n", window -> ibs[INDEX(i, j, N)].ibs0, window -> ibs[INDEX(i, j, N)].ibs1, window -> ibs[INDEX(i, j, N)].ibs2);
                    }
                }
            } else {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j <= i; j++) {
                        if (asdToIbs)
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
        fprintf(output, "\n");
    }
}

// Defines all the command line options supplied
//  by user and parameters to run the LODESTAR analysis.
typedef struct {
    // The input VCF file.
    char* input_filename;
    // The parser object used to read the VCF file.
    VCFLocusParser_t* parser;
    // The basename for output files.
    char* output_basename;
    // Pointer to the file summarizing sliding window.
    FILE* windows_info;
    // Pointer to the file holding asd matrices and points of
    //  samples in reduced space.
    FILE* windows_coords;
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
    char* target_filename;
    // Pointer to file of user specified points to perform Procrustes
    //  analysis with.
    FILE* target_file;
    // Pvalue threshold used to print points.
    double pthresh;
    // Number of permutations to execute in permutation test.
    int num_perms;
    // t-statistic threshold used to print points.
    double tthresh;
    // Regions to include/exclude while parsing VCF file.
    char* regions;
    // Minor allele frequency threshhold.
    double maf;
    // Missing allele frequency threshold.
    double afMissing;
    // Regions to save window IBS/ASD matrix.
    char* save_regions_str;
    RegionSet_t* save_regions;
    // Covnert ASD values to IBS counts.
    bool asdToIbs;
    // Regions to print coodinates of samples in dimension-k.
    char* print_regions_str;
    RegionSet_t* print_regions;
    // Output long format instead of lower triangle for pairwise calculations.
    bool long_output;
    // Output file in JSON format instead of text file.
    bool json_output;
} LodestarConfiguration_t;

// Make sure user defined arguments are valid.
// Accepts:
//  LodestarConfiguration_t* lodestar_config -> The configured parameters.
// Returns: int, 0 if all parameters are valid. 1 if user supplied an invalid value.
int check_configuration(LodestarConfiguration_t* lodestar_config) {
    // Check maf.
    if (lodestar_config -> maf < 0 || lodestar_config -> maf > 1) { fprintf(stderr, "--maf must be in [0, 1].\n"); return 1;}
    // Check afMissing.
    if (lodestar_config -> afMissing < 0 || lodestar_config -> afMissing > 1) { fprintf(stderr, "--afMissing must be in [0, 1].\n"); return 1;}
    // Check that input VCF file exists and region is valid if supplied.
    kstring_t* r = (kstring_t*) calloc(1, sizeof(kstring_t));
    bool takeComplement; 
    if (lodestar_config -> regions != NULL) {
        takeComplement = lodestar_config -> regions[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestar_config -> regions + 1, r);
        else 
            ks_overwrite(lodestar_config -> regions, r);
        if ((lodestar_config -> parser = init_vcf_locus_parser(lodestar_config -> input_filename, r, takeComplement, lodestar_config -> maf, lodestar_config -> afMissing, true)) == NULL) {
            free(ks_str(r)); free(r);
            fprintf(stderr, "-i file does not exist or --regions invalid.\n");
            return 1;
        }
    } else {
        if ((lodestar_config -> parser = init_vcf_locus_parser(lodestar_config -> input_filename, NULL, false, lodestar_config -> maf, lodestar_config -> afMissing, true)) == NULL) {
            free(r);
            fprintf(stderr, "-i file does not exist.\n");
            return 1;
        }
    }
    // Check regions to save is valid.
    if (lodestar_config -> save_regions_str != NULL) {
        takeComplement = lodestar_config -> save_regions_str[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestar_config -> save_regions_str + 1, r);
        else 
            ks_overwrite(lodestar_config -> save_regions_str, r);
        if ((lodestar_config -> save_regions = init_region_set(r, takeComplement)) == NULL) {
            free(ks_str(r)); free(r);
            fprintf(stderr, "--save region is invalid.\n");
            return 1;
        }
    }
    // Check region to print is valid.
    if (lodestar_config -> print_regions_str != NULL) {
        takeComplement = lodestar_config -> print_regions_str[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestar_config -> print_regions_str + 1, r);
        else 
            ks_overwrite(lodestar_config -> print_regions_str, r);
        if ((lodestar_config -> print_regions = init_region_set(r, takeComplement)) == NULL) {
            free(ks_str(r)); free(r);
            fprintf(stderr, "--printCoords region is invalid.\n");
            return 1;
        }
    }
    if (ks_str(r) != NULL)
        free(ks_str(r)); 
    free(r);
    // Check that output basename was given.
    if (lodestar_config -> output_basename == NULL) { fprintf(stderr, "-o was not given.\n"); return 1;}
    // Check haplotype size is valid.
    if (lodestar_config -> HAP_SIZE <= 0) { fprintf(stderr, "-h must be INT > 0.\n"); return 1;}
    // Check window size is valid.
    if (!lodestar_config -> global && lodestar_config -> WINDOW_SIZE <= 0) { fprintf(stderr, "-w must be INT > 0 supplied by the user.\n"); return 1;}
    // Check step size is valid.
    if (!lodestar_config -> global && (lodestar_config -> STEP_SIZE <= 0 || lodestar_config -> STEP_SIZE > lodestar_config -> WINDOW_SIZE)) { fprintf(stderr, "-s must be 0 < INT <= WINDOW_SIZE supplied by the user.\n"); return 1;}
    // Check k.
    if (lodestar_config -> k <= 0 || lodestar_config -> k >= lodestar_config -> parser -> numSamples) { fprintf(stderr, "-k must be 1 <= INT < N.\n"); return 1;}
    // Check number of threads.
    if (lodestar_config -> threads <= 0) { fprintf(stderr, "--threads must be INT > 0.\n"); return 1;}
    // Check pthresh.
    if (lodestar_config -> pthresh < 0 || lodestar_config -> pthresh > 1) { fprintf(stderr, "--pthresh must be in (0, 1].\n"); return 1;}
    // Check number of permutation.
    if (lodestar_config -> num_perms <= 0) { fprintf(stderr, "--perms must be INT > 0.\n"); return 1;}
    // Check tthresh.
    if (lodestar_config -> tthresh <= 0 || lodestar_config -> tthresh > 1) { fprintf(stderr, "--tthresh must be in (0, 1].\n"); return 1;}
    // Check target file if exists.
    if (lodestar_config -> target_filename != NULL) {
        lodestar_config -> target_file = fopen(lodestar_config -> target_filename, "r");
        if (lodestar_config -> target_file == NULL) {
            fprintf(stderr, "--target %s does not exist.\n", lodestar_config -> target_filename);
            return 1;
        }
        int numLines = 1, dim = 1;
        char c;
        while ((c = getc(lodestar_config -> target_file)) != EOF) {
            if (numLines == 0 && c == ',')
                dim++;
            if (c == '\n')
                numLines++;
        }
        fclose(lodestar_config -> target_file);
        // Check the number of samples in target file.
        if (numLines != lodestar_config -> parser -> numSamples) {
            fprintf(stderr, "Number of samples in --target %s does not match that of the input file.\n", lodestar_config -> target_filename);
            return 1;
        }
        // Check dimension of target file.
        if (dim != lodestar_config -> k) {
            fprintf(stderr, "Dimension of --target %s does not match that of k.\n", lodestar_config -> target_filename);
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
void print_configuration(FILE* output, LodestarConfiguration_t lodestar_config) {
    fprintf(output, "LODESTAR Configuration\n");
    fprintf(output, "----------------------\n");
    fprintf(output, "Input file: %s\n", lodestar_config.input_filename);
    fprintf(output, "Output basename: %s\n", lodestar_config.output_basename);
    fprintf(output, "Haplotype size: %d\n", lodestar_config.HAP_SIZE);
    fprintf(output, "Window size: %d\n", lodestar_config.WINDOW_SIZE);
    fprintf(output, "Step size: %d\n", lodestar_config.STEP_SIZE);
    fprintf(output, "K: %d\n", lodestar_config.k);
    fprintf(output, "Use similarity: %s\n", PRINT_BOOL(lodestar_config.similarity));
    fprintf(output, "Calculate genome-wide only: %s\n", PRINT_BOOL(lodestar_config.global));
    fprintf(output, "Number of threads used: %d\n", lodestar_config.threads);
    fprintf(output, "P-value threshold (Disabled if 0): %lf\n", lodestar_config.pthresh);
    fprintf(output, "Number of permutations (Disabled p-value threshold is 0): %d\n", lodestar_config.num_perms);
    fprintf(output, "Procrustes statistic threshold: %lf\n", lodestar_config.tthresh);
    if (lodestar_config.regions != NULL)
        fprintf(output, "Include/exclude records from input: %s\n", lodestar_config.regions);
    fprintf(output, "Minor allele frequency threshold: %lf\n", lodestar_config.maf);
    fprintf(output, "Missing allele frequency threshold: %lf\n", lodestar_config.afMissing);
    fprintf(output, "Convert ASD values to IBS counts: %s\n", PRINT_BOOL(lodestar_config.asdToIbs));
    if (lodestar_config.save_regions_str != NULL)
        fprintf(output, "Save IBS/ASD values for windows overlapping/not overlapping: %s\n", lodestar_config.save_regions_str);
    if (lodestar_config.print_regions_str != NULL)
        fprintf(output, "Print coordinates of windows that are/are not overlapping: %s\n", lodestar_config.print_regions_str);
    if (lodestar_config.target_filename != NULL)
        fprintf(output, "File of coordinates to perform Procrustes analysis against: %s\n", lodestar_config.target_filename);
    fprintf(output, "Use long-format output: %s\n", PRINT_BOOL(lodestar_config.long_output));
    fprintf(output, "Save as JSON file: %s\n", PRINT_BOOL(lodestar_config.json_output));
}

// Destroy all dynamically allocated memory associated with lodestar_config.
// Accepts:
//  LodestarConfiguration_t lodestar_config -> The configuration to destroy.
// Returns: void.
void destroy_lodestar_configuration(LodestarConfiguration_t lodestar_config) {
    if (lodestar_config.parser != NULL)
        destroy_vcf_locus_parser(lodestar_config.parser);
    if (lodestar_config.windows_info != NULL)
        fclose(lodestar_config.windows_info);
    if (lodestar_config.windows_coords != NULL)
        fclose(lodestar_config.windows_coords);
    if (lodestar_config.target_file != NULL)
        fclose(lodestar_config.target_file);
    if (lodestar_config.print_regions != NULL)
        destroy_region_set(lodestar_config.print_regions);
    if (lodestar_config.save_regions != NULL)
        destroy_region_set(lodestar_config.save_regions);
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
    fprintf(stderr, "Usage: lodestar [options] -w <WINDOW_SIZE> -s <STEP_SIZE> -i <input.vcf.gz> -o <output_basename>\n");
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
    fprintf(stderr, "  --regions [^]REGS       Include/exclude records from VCF defined by REGS.\n");
    fprintf(stderr, "                               Default NULL.\n");
    fprintf(stderr, "   --printCoords [^]REGS   Print coordinates that are overlapping/non-overlapping with REGS.\n");
    fprintf(stderr, "                               Default NULL.\n");
    fprintf(stderr, "   --save [^]REGS          Save IBS/ASD values for windows that are overlapping/non-overlapping with REGS.\n");
    fprintf(stderr, "                               Default NULL. See --asdToIbs\n");
    fprintf(stderr, "   --asdToIbs              Convert IBS values to ASD values in output.\n");
    fprintf(stderr, "                               Default false.\n");
    fprintf(stderr, "   --maf DOUBLE            Drops VCF records with a MAF less than threshold.\n");
    fprintf(stderr, "                               Default 0.\n");
    fprintf(stderr, "   --afMissing DOUBLE      Drops VCF records with fraction of missing genotypes greater than or equal to threshold.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   --ibs [^]REGS           Saves IBS calculations for overlapping windows included/excluded defined by REGS.\n");
    fprintf(stderr, "                               Not used when --global is set.\n");
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
    //  Also, check if any options are suplpied that are not defined. If user supplies the help
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
    LodestarConfiguration_t lodestar_config;
    lodestar_config.input_filename = NULL;
    lodestar_config.parser = NULL;
    lodestar_config.output_basename = NULL;
    lodestar_config.windows_info = NULL;
    lodestar_config.windows_coords = NULL;
    lodestar_config.HAP_SIZE = 1;
    lodestar_config.WINDOW_SIZE = -1;
    lodestar_config.STEP_SIZE = -1;
    lodestar_config.k = 2;
    lodestar_config.threads = 1;
    lodestar_config.similarity = false;
    lodestar_config.global = false;
    lodestar_config.target_filename = NULL;
    lodestar_config.target_file = NULL;
    lodestar_config.pthresh = 0;
    lodestar_config.num_perms = 10000;
    lodestar_config.tthresh = 0.95;
    lodestar_config.regions = NULL;
    lodestar_config.maf = 0;
    lodestar_config.afMissing = 1;
    lodestar_config.save_regions_str = NULL;
    lodestar_config.save_regions = NULL;
    lodestar_config.asdToIbs = false;
    lodestar_config.print_regions_str = NULL;
    lodestar_config.print_regions = NULL;
    lodestar_config.long_output = false;
    lodestar_config.json_output = false;

    // Parse command line arguments.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'i': lodestar_config.input_filename = options.arg; break;
            case 'o': lodestar_config.output_basename = options.arg; break;
            case 'h': lodestar_config.HAP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'w': lodestar_config.WINDOW_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 's': lodestar_config.STEP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'k': lodestar_config.k = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 302: lodestar_config.similarity = true; break;
            case 303: lodestar_config.global = true; break;
            case 304: lodestar_config.threads = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 305: lodestar_config.pthresh = strtod(options.arg, (char**) NULL); break;
            case 306: lodestar_config.num_perms = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 307: lodestar_config.tthresh = strtod(options.arg, (char**) NULL); break;
            case 308: lodestar_config.regions = options.arg; break;
            case 309: lodestar_config.maf = strtod(options.arg, (char**) NULL); break;
            case 310: lodestar_config.afMissing = strtod(options.arg, (char**) NULL); break;
            case 311: lodestar_config.save_regions_str = options.arg; break;
            case 312: lodestar_config.asdToIbs = true; break;
            case 313: lodestar_config.long_output = true; break;
            case 314: lodestar_config.json_output = true; break;
            case 315: lodestar_config.target_filename = options.arg; break;
            case 316: lodestar_config.print_regions_str = options.arg; break;
        }
	}
    
    // Check configuration. If invalid argument, exit program.
    if (check_configuration(&lodestar_config) != 0) {
        fprintf(stderr, "Exiting!\n");
        return 1;
    }

    // If this point is reached, the configuration is valid, and
    //  we can execute our analysis.

    // Setup the logfile.
    kstring_t* outputBasename = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs(lodestar_config.output_basename, outputBasename);
    kputs(".log", outputBasename);
    INIT_LOG(ks_str(outputBasename));
    // Print configuration to log file.
    print_configuration(log -> file, lodestar_config);
    // Print sample names in log file.
    fprintf(log -> file, "\nSample Names:\n");
    for (int i = 0; i < lodestar_config.parser -> numSamples; i++)
        fprintf(log -> file, "%s\n", lodestar_config.parser -> sampleNames[i].s);
    fprintf(log -> file, "\n");
    printf("Logging progress in %s\n", ks_str(outputBasename));



    // Close all files and free all memory used in analysis.
    printf("Done!\n");
    CLOSE_LOG();
    free(ks_str(outputBasename)); free(outputBasename);
    destroy_lodestar_configuration(lodestar_config);
    return 0;

}