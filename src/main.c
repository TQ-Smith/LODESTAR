
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

typedef struct {
    char* input_filename;
    VCFLocusParser_t* parser;
    char* output_basename;
    FILE* windows_info;
    FILE* windows_coords;
    int HAP_SIZE;
    int WINDOW_SIZE;
    int STEP_SIZE;
    int k;
    int threads;
    bool similarity;
    bool global;
    char* target_filename;
    FILE* target_file;
    double pthresh;
    int num_perms;
    double tthresh;
    char* regions;
    double maf;
    double afMissing;
    char* ibs_regions_str;
    RegionSet_t* ibs_regions;
    char* asd_regions_str;
    RegionSet_t* asd_regions;
    char* print_regions_str;
    RegionSet_t* print_regions;
    bool long_output;
    bool json_output;
} LodestarConfiguration;

int check_configuration(LodestarConfiguration* lodestar_config) {
    
    if (lodestar_config -> maf < 0 || lodestar_config -> maf > 1) {
        fprintf(stderr, "MAF must be in [0, 1].\n");
        return 1;
    }

    if (lodestar_config -> afMissing < 0 || lodestar_config -> afMissing > 1) {
        fprintf(stderr, "Missing genotype threshold must be in [0, 1].\n");
        return 1;
    }
    
    kstring_t* r = (kstring_t*) calloc(1, sizeof(kstring_t));
    bool takeComplement; 
    if (lodestar_config -> regions != NULL) {
        takeComplement = lodestar_config -> regions[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestar_config -> regions + 1, r);
        else 
            ks_overwrite(lodestar_config -> regions, r);
    }
    if ((lodestar_config -> parser = init_vcf_locus_parser(lodestar_config -> input_filename, r, takeComplement, lodestar_config -> maf, lodestar_config -> afMissing, true)) == NULL) {
        free(ks_str(r)); free(r);
        fprintf(stderr, "Invalid input file name or region.\n");
        return 1;
    }

    if (lodestar_config -> ibs_regions_str != NULL) {
        takeComplement = lodestar_config -> ibs_regions_str[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestar_config -> ibs_regions_str + 1, r);
        else 
            ks_overwrite(lodestar_config -> ibs_regions_str, r);
    }
    if ((lodestar_config -> ibs_regions = init_region_set(r, takeComplement)) == NULL) {
        free(ks_str(r)); free(r);
        fprintf(stderr, "Invlaid region given for IBS option.\n");
        return 1;
    }

    if (lodestar_config -> asd_regions_str != NULL) {
        takeComplement = lodestar_config -> asd_regions_str[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestar_config -> asd_regions_str + 1, r);
        else 
            ks_overwrite(lodestar_config -> asd_regions_str, r);
    }
    if ((lodestar_config -> asd_regions = init_region_set(r, takeComplement)) == NULL) {
        free(ks_str(r)); free(r);
        fprintf(stderr, "Invlaid region given for ASD option.\n");
        return 1;
    }

    if (lodestar_config -> print_regions_str != NULL) {
        takeComplement = lodestar_config -> print_regions_str[0] == '^' ? true : false;
        if (takeComplement)
            ks_overwrite(lodestar_config -> print_regions_str + 1, r);
        else 
            ks_overwrite(lodestar_config -> print_regions_str, r);
    }
    if ((lodestar_config -> print_regions = init_region_set(r, takeComplement)) == NULL) {
        free(ks_str(r)); free(r);
        fprintf(stderr, "Invlaid region given for printCoords option.\n");
        return 1;
    }
    free(ks_str(r)); free(r);

    if (lodestar_config -> output_basename == NULL) {
        fprintf(stderr, "No output basename given.\n");
        return 1;
    }
    
    if (lodestar_config -> HAP_SIZE <= 0) {
        fprintf(stderr, "Haplotype size must be a positive integer.\n");
        return 1;
    }

    if (lodestar_config -> WINDOW_SIZE <= 0) {
        fprintf(stderr, "Window size must be a positive integer supplied by the user.\n");
        return 1;
    }

    if (lodestar_config -> STEP_SIZE <= 0 || lodestar_config -> STEP_SIZE > lodestar_config -> WINDOW_SIZE) {
        fprintf(stderr, "Step size must be a positive integer less than or equal to the window size supplied by the user.\n");
        return 1;
    }

    if (lodestar_config -> k <= 0 || lodestar_config -> k >= lodestar_config -> parser -> numSamples) {
        fprintf(stderr, "k must be an integer between 1 and (number of samples - 1).\n");
        return 1;
    }

    if (lodestar_config -> threads <= 0) {
        fprintf(stderr, "threads must be a positive integer greater than 0.\n");
        return 1;
    }

    if (lodestar_config -> pthresh <= 0 || lodestar_config -> pthresh > 1) {
        fprintf(stderr, "pthresh must be a real number greater than 0 and less than or equal to 1.\n");
        return 1;
    }

    if (lodestar_config -> num_perms <= 0) {
        fprintf(stderr, "perms must be a positive integer greater than 0.\n");
        return 1;
    }

    if (lodestar_config -> tthresh <= 0 || lodestar_config -> tthresh > 1) {
        fprintf(stderr, "tthresh must be a real number greater than 0 and less than or equal to 1.\n");
        return 1;
    }

    if (lodestar_config -> target_filename != NULL) {
        lodestar_config -> target_file = fopen(lodestar_config -> target_filename, "r");
        if (lodestar_config -> target_file == NULL) {
            fprintf(stderr, "%s does not exist.\n", lodestar_config -> target_filename);
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
        if (numLines != lodestar_config -> parser -> numSamples) {
            fprintf(stderr, "Number of samples in coordinate file %s does not match that of the input file.\n", lodestar_config -> target_filename);
            return 1;
        }
        if (dim != lodestar_config -> k) {
            fprintf(stderr, "Dimension of coordinate file %s does not match that of k.\n", lodestar_config -> target_filename);
            return 1;
        }
    }

    return 0;
}

#define PRINT_BOOL(X) (X ? "true" : "false")

void print_configuration(FILE* output, LodestarConfiguration lodestar_config) {
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
    if (lodestar_config.regions != NULL && lodestar_config.regions[0] == '^')
        fprintf(output, "Exclude records from input: %s\n", lodestar_config.regions);
    if (lodestar_config.regions != NULL && lodestar_config.regions[0] != '^')
        fprintf(output, "Include records from input: %s\n", lodestar_config.regions);
    fprintf(output, "Minor allele frequency threshold: %lf\n", lodestar_config.maf);
    fprintf(output, "Missing allele frequency threshold: %lf\n", lodestar_config.afMissing);
    if (lodestar_config.ibs_regions_str != NULL && lodestar_config.ibs_regions_str[0] == '^')
        fprintf(output, "Exclude records for saving IBS values: %s\n", lodestar_config.ibs_regions_str);
    if (lodestar_config.ibs_regions_str != NULL && lodestar_config.ibs_regions_str[0] != '^')
        fprintf(output, "Include records for saving IBS values: %s\n", lodestar_config.ibs_regions_str);
    if (lodestar_config.asd_regions_str != NULL && lodestar_config.asd_regions_str[0] == '^')
        fprintf(output, "Exclude records for saving ASD values: %s\n", lodestar_config.asd_regions_str);
    if (lodestar_config.asd_regions_str != NULL && lodestar_config.asd_regions_str[0] != '^')
        fprintf(output, "Include records for saving ASD values: %s\n", lodestar_config.asd_regions_str);
    if (lodestar_config.print_regions_str != NULL && lodestar_config.print_regions_str[0] == '^')
        fprintf(output, "Print coordinates of windows that are not overlapping: %s\n", lodestar_config.print_regions_str);
    if (lodestar_config.print_regions_str != NULL && lodestar_config.print_regions_str[0] != '^')
        fprintf(output, "Print coordinates of windows that are overlapping: %s\n", lodestar_config.print_regions_str);
    if (lodestar_config.target_filename != NULL)
        fprintf(output, "File of coordinates to perform Procrustes analysis against: %s\n", lodestar_config.target_filename);
    fprintf(output, "Use long-format output: %s\n", PRINT_BOOL(lodestar_config.long_output));
    fprintf(output, "Save as JSON file: %s\n", PRINT_BOOL(lodestar_config.json_output));
}

void destroy_lodestar_configuration(LodestarConfiguration lodestar_config) {
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
    if (lodestar_config.ibs_regions != NULL)
        destroy_region_set(lodestar_config.ibs_regions);
    if (lodestar_config.asd_regions != NULL)
        destroy_region_set(lodestar_config.asd_regions);
}

void print_help() {
    fprintf(stderr, "\n");
    fprintf(stderr, "LODESTAR v1.0 May 2024\n");
    fprintf(stderr, "----------------------\n\n");
    fprintf(stderr, "Written by T. Quinn Smith\n");
    fprintf(stderr, "Principal Investigator: Zachary A. Szpiech\n");
    fprintf(stderr, "The Pennsylvania State University\n\n");
    fprintf(stderr, "Usage: lodestar [options] -i <input.vcf.gz> -o <output_basename>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   --help                      Prints help menu and exits.\n");
    fprintf(stderr, "   --version                   Prints version number and exits.\n");
    fprintf(stderr, "   -i STR                      Path to input vcf file.\n");
    fprintf(stderr, "   -o STR                      Output base names for files.\n");
    fprintf(stderr, "   -h INT                      Number of loci in a haplotype.\n");
    fprintf(stderr, "                                   Default 1. Not used when --global is set.\n");
    fprintf(stderr, "   -w INT                      Number of haplotypes in a window.\n");
    fprintf(stderr, "                                   Must be set by user. Not used when --global is set.\n");
    fprintf(stderr, "   -s INT                      Number of haplotypes to increment the sliding window.\n");
    fprintf(stderr, "                                   Must be set by user. Not used when --global is set.\n");
    fprintf(stderr, "   -k INT                      Dimension to project samples into. Must be less than number of samples.\n");
    fprintf(stderr, "                                   Default 2. Must be less than number of samples in VCF.\n");
    fprintf(stderr, "   --threads INT               Number of threads to use in computation.\n");
    fprintf(stderr, "                                   Default 1.\n");
    fprintf(stderr, "   --similarity                Compute similarity between sets of points instead of dissimilarity.\n");
    fprintf(stderr, "   --global                    Compute only the global set of points.\n");
    fprintf(stderr, "                                   Takes presedence over windowing parameters.\n");
    fprintf(stderr, "   --target STR                A n-by-k csv file containing user defined coordinates to perform Procrustes analysis.\n");
    fprintf(stderr, "   --pthresh DOUBLE            Print coordinates of all windows less than or equal to threshold.\n");
    fprintf(stderr, "                                   Window must also satisfy tthresh. Default 0.\n");
    fprintf(stderr, "   --perms INT                 The number of permutations to execute.\n");
    fprintf(stderr, "                                   Default 10000. Permutation test does not execute if pthresh is 0.\n");
    fprintf(stderr, "   --tthresh DOUBLE            Print coordinates of all windows greater than or equal to threshold.\n");
    fprintf(stderr, "                                   Window must also satisfy pthresh. Default 0.95.\n");
    fprintf(stderr, "   --printCoords [^]REGIONS    Print coordinates that are overlapping/non-overlapping with REGIONS regardless of pthresh and tthresh.\n");
    fprintf(stderr, "   --regions [^]REGIONS        Include/exclude records from VCF defined by REGIONS.\n");
    fprintf(stderr, "   --maf DOUBLE                Drops VCF records with a MAF less than threshold.\n");
    fprintf(stderr, "                                   Default 0.\n");
    fprintf(stderr, "   --afMissing DOUBLE          Drops VCF records with fraction of missing genotypes greater than or equal to threshold.\n");
    fprintf(stderr, "                                   Default 1.\n");
    fprintf(stderr, "   --ibs [^]REGIONS            Saves IBS calculations for overlapping windows included/excluded defined by REGIONS.\n");
    fprintf(stderr, "                                   Not used when --global is set.\n");
    fprintf(stderr, "   --asd [^]REGIONS            Saves ASD calculations for overlapping windows included/excluded defined by REGIONS.\n");
    fprintf(stderr, "                                   Not used when --global is set.\n");
    fprintf(stderr, "   --long                      Prints calculations in long format instead of matrix form.\n");
    fprintf(stderr, "   --json                      Prints window information in JSON format instead of TXT.\n");
    fprintf(stderr, "Types:\n");
    fprintf(stderr, "   STR                     A string.\n");
    fprintf(stderr, "   INT                     A non-negative integer.\n");
    fprintf(stderr, "   DOUBLE                  A real number between 0 and 1, inclusive.\n");
    fprintf(stderr, "   REGIONS                 REGION,REGIONS | REGION\n");
    fprintf(stderr, "   REGION                  STR | STR:INT | STR:-INT | STR:INT- | STR:INT-INT\n");
    fprintf(stderr, "\n");
}

static ko_longopt_t long_options[] = {
    {"help",            ko_no_argument,         300},
    {"version",         ko_no_argument,         'v'},
    {"similarity",      ko_no_argument,         302},
    {"global",          ko_no_argument,         303},
    {"threads",         ko_required_argument,   304},
    {"pthresh",         ko_required_argument,   305},
    {"perms",           ko_required_argument,   306},
    {"tthresh",         ko_required_argument,   307},
    {"[^]regions",      ko_required_argument,   308},
    {"maf",             ko_required_argument,   309},
    {"afMissing",       ko_required_argument,   310},
    {"[^]ibs",          ko_required_argument,   311},
    {"[^]asd",          ko_required_argument,   312},
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

    const char *opt_str = "i:o:h:w:s:k:v";
    ketopt_t options = KETOPT_INIT;
    int c;

    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option %s is missing an argument! Exiting ...\n", argv[options.i - 1]); return 1;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return 1;
            case 300: print_help(); return 0;
            case 'v': fprintf(stderr, "Version 1.0 May 2024.\n"); return 0;
        }
	}
	options = KETOPT_INIT;

    LodestarConfiguration lodestar_config;

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
    lodestar_config.ibs_regions_str = NULL;
    lodestar_config.ibs_regions = NULL;
    lodestar_config.asd_regions_str = NULL;
    lodestar_config.asd_regions = NULL;
    lodestar_config.print_regions_str = NULL;
    lodestar_config.print_regions = NULL;
    lodestar_config.long_output = false;
    lodestar_config.json_output = false;

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
            case 311: lodestar_config.ibs_regions_str = options.arg; break;
            case 312: lodestar_config.asd_regions_str = options.arg; break;
            case 313: lodestar_config.long_output = true; break;
            case 314: lodestar_config.json_output = true; break;
            case 315: lodestar_config.target_filename = options.arg; break;
            case 316: lodestar_config.print_regions_str = options.arg; break;
        }
	}
    
    if (check_configuration(&lodestar_config) != 0) {
        fprintf(stderr, "Exiting!\n");
        return 1;
    }

    return 0;

}