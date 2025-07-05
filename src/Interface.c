
// File: Interface.c
// Date: 18 February 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines LODESTAR configuration from command line interface and
//              output formatting of LODESTAR results.

#include "Interface.h"
#include "ketopt.h"
#include <unistd.h> // Check if file exits.

// Make sure user defined arguments are valid.
// Accepts:
//  LodestarConfig_t* lodestarConfig -> The configured parameters.
// Returns: int, 0 if all parameters are valid. -1 if user supplied an invalid value.
int check_configuration(LodestarConfig_t* lodestarConfig) {
    // Check maf.
    if (lodestarConfig -> maf < 0 || lodestarConfig -> maf > 1) { 
        fprintf(stderr, "--maf %lf must be in [0, 1].\n", lodestarConfig -> maf); 
        return -1;
    }
    // Check afMissing.
    if (lodestarConfig -> afMissing < 0 || lodestarConfig -> afMissing > 1) { 
        fprintf(stderr, "--afMissing %lf must be in [0, 1].\n", lodestarConfig -> afMissing); 
        return -1;
    }
    if (access(lodestarConfig -> inputFileName, F_OK) != 0) {
        fprintf(stderr, "Input file %s does not exist.\n", lodestarConfig -> inputFileName);
        return -1;
    }
    // Check target file if exists.
    if (lodestar_config -> targetFileName != NULL && access(lodestar_config -> targetFileName, F_OK) != 0) {
        fprintf(stderr, "-y %s does not exist.\n", lodestar_config -> targetFileName); 
        return -1;
    }
    // Check haplotype and block size are valid.
    if (lodestarConfig -> BLOCK_SIZE <= 0) { 
        fprintf(stderr, "-b %d must be INT > 0.\n", lodestarConfig -> BLOCK_SIZE); 
        return -1;
    }
    if (lodestarConfig -> HAP_SIZE <= 0 || lodestarConfig -> HAP_SIZE > lodestarConfig -> BLOCK_SIZE) { 
        fprintf(stderr, "-h %d must be 0 < INT < BLOCK_SIZE.\n", lodestarConfig -> HAP_SIZE); 
        return -1;
    }
    // Check number of threads.
    if (lodestarConfig -> threads <= 0) { 
        fprintf(stderr, "--threads %d must be INT > 0.\n", lodestarConfig -> threads); 
        return -1;
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
    fprintf(stderr, "Usage: lodestar [options] <input>.vcf.gz\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   -o STR                  Output file basename.\n");
    fprintf(stderr, "                               Default <input>.\n");
    fprintf(stderr, "   -h INT                  Number of loci in a haplotype.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   -b INT                  The block size in base pairs.\n");
    fprintf(stderr, "                               Default 2MB.\n");
    fprintf(stderr, "   -k INT                  Dimension to project samples into. Must be less than number of samples.\n");
    fprintf(stderr, "                               Default 2. Must be less than number of samples in VCF.\n");
    fprintf(stderr, "   -y file.tsv             A n-by-k tsv file containing user defined coordinates.\n");
    fprintf(stderr, "                               The set of points to use in Procrustes analysis\n");
    fprintf(stderr, "   -t INT                  Number of threads to use in computation.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   --maf DOUBLE            Drops biallelic VCF records with a MAF less than threshold.\n");
    fprintf(stderr, "                               Default 0.05.\n");
    fprintf(stderr, "   --afMissing DOUBLE      Drops VCF records with fraction of genotypes missing greater than threshold.\n");
    fprintf(stderr, "                               Default 0.1\n");
    fprintf(stderr, "\n");
}

// Long options used for LODESTAR.
static ko_longopt_t long_options[] = {
    {"threads",         ko_required_argument,   't'},
    {"maf",             ko_required_argument,   309},
    {"afMissing",       ko_required_argument,   310},
    {"target",          ko_required_argument,   'y'},
    {"dimension",       ko_required_argument,   'k'},
    {"haplotypeSize",   ko_required_argument,   'h'},
    {"blockSize",       ko_required_argument,   'b'},
    {0, 0, 0}
};

LodestarConfig_t* init_lodestar_config(int argc, char *argv[]) {

    // Print help menu when no options/arguments were given.
    if (argc == 1) {
        print_help();
        return NULL;
    } 

    // Parse user options.
    const char *opt_str = "h:b:k:t:y:";
    ketopt_t options = KETOPT_INIT;
    int c;

    // Pass through options. Check for options requiring an argument that were not given one.
    //  Also, check if any options are supplied that are not defined. If user supplies the help
    //  argument, print help menu and exit, likewise for version.
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return NULL;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return NULL;
        }
	}

    // Set defaults for LODESTAR configuration.
    LodestarConfig_t* lodestarConfig = calloc(1, sizeof(LodestarConfig_t));
    lodestarConfig -> inputFileName = NULL;
    lodestarConfig -> outputBasename = NULL;
    lodestarConfig -> HAP_SIZE = 1;
    lodestarConfig -> BLOCK_SIZE = 2000000;
    lodestarConfig -> k = 2;
    lodestarConfig -> threads = 1;
    lodestarConfig -> targetFileName = NULL;
    lodestarConfig -> maf = 0.05;
    lodestarConfig -> afMissing = 0.1;

    // Parse command line arguments.
    options = KETOPT_INIT;
    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case 'o': lodestarConfig -> outputBasename = strdup(options.arg); break;
            case 'h': lodestarConfig -> HAP_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'b': lodestarConfig -> BLOCK_SIZE = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'k': lodestarConfig -> k = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 't': lodestarConfig -> threads = (int) strtol(options.arg, (char**) NULL, 10); break;
            case 'y': lodestarConfig -> targetFileName = strdup(options.arg); break;
            case 309: lodestarConfig -> maf = strtod(options.arg, (char**) NULL); break;
            case 310: lodestarConfig -> afMissing = strtod(options.arg, (char**) NULL); break;
        }
	}

    // If the input file was not give, exit.
    if (argc - options.ind != 1) {
        fprintf(stderr, "No input file given. Exiting!\n");
        destroy_lodestar_config(config);
        return NULL;
    }

    // Set the input filename.
    lodestarConfig -> inputFileName = strdup(argv[options.ind]);

    // Check to that user input is all valid.
    if (check_configuration(lodestarConfig) != 0)
        return NULL;

    // If output base name was not specific, then use the input file basename.
    if (config -> outBaseName == NULL) {
        int endPos = (int) ((long) strstr(config -> inputFileName, ".vcf") - (long) config -> inputFileName) + 1;
        config -> outBaseName = strndup(config -> inputFileName, endPos - 1);
    }

    // Copy command to echo in output files.
    kstring_t* cmd = calloc(1, sizeof(kstring_t));
    for (int i = 0; i < argc; i++)
        kvsprintf(cmd, "%s", argv[i]);
    lodestarConfig -> cmd = strdup(cmd -> s);
    free(cmd -> s); free(cmd);

    return lodestarConfig;

}

void destroy_lodestar_config(LodestarConfig_t* lodestarConfig) {
    if (lodestarConfig == NULL)
        return NULL;
    if (lodestarConfig -> targetFileName != NULL)
        free(lodestarConfig -> targetFileName);
    free(lodestarConfig -> inputFileName);
    free(lodestarConfig -> outputBasename);
    free(lodestarConfig);
}