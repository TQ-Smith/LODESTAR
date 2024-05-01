
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include "ProcrustesAnalysis.h"

#include "../lib/ketopt.h"

typedef struct {
    VCFLocusParser* parser;
    char* output_basename;
    FILE* windowsInfo;
    FILE* windowsCoords;
    int HAP_SIZE;
    int WINDOW_SIZE;
    int STEP_SIZE;
    int k;
    int threads;
    bool similarity;
    bool global;
    char* coords_file_name;
    FILE* coords_file;
    double pthresh;
    int num_perms;
    double tthresh;
    char* regions;
    double maf;
    double afMissing;
    char* ibs_regions_str;
    RegionFilter* ibs_regions;
    char* asd_regions_str;
    RegionFilter* asd_regions;
    bool long_output;
    bool json_output;
} LodestarConfiguration;

void print_configuration(FILE* output, LodestarConfiguration lodestar_config) {
    fprintf(output, "Parameters used to run LODESTAR:\n");
    fprintf(output, "--------------------------------\n");
    if (lodestar_config.parser != NULL)
        fprintf(output, "Input File: %s\n", ks_str(lodestar_config.parser -> fileName));
    if (lodestar_config.output_basename != NULL) {
        if (!lodestar_config.global)
            fprintf(output, "Window Statistics File: %s%s.txt\n", lodestar_config.output_basename, "Statistics");
        fprintf(output, "Window Coordinates File: %s%s.%s\n", lodestar_config.output_basename, "Coordinates", lodestar_config.json_output ? "json" : "txt");
    }
    
}

void destroy_lodestar_configuration(LodestarConfiguration lodestar_config) {

}

void print_help() {
    fprintf(stderr, "LODESTAR v1.0 May 2024\n");
    fprintf(stderr, "----------------------\n\n");
    fprintf(stderr, "Written by T. Quinn Smith\n");
    fprintf(stderr, "Principal Investigator: Zachary A. Szpiech\n");
    fprintf(stderr, "The Pennsylvania State University\n\n");
    fprintf(stderr, "Usage: lodestar [options] -i <input.vcf.gz> -o <output_basename>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   --help                  Prints help menu and exits.\n");
    fprintf(stderr, "   --version               Prints version number and exits.\n");
    fprintf(stderr, "   -i STR                  Path to input vcf file.\n");
    fprintf(stderr, "   -o STR                  Output base names for files.\n");
    fprintf(stderr, "   -h INT                  Number of loci in a haplotype.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   -w INT                  Number of haplotypes in a window.\n");
    fprintf(stderr, "                               Default 100.\n");
    fprintf(stderr, "   -s INT                  Number of haplotypes to increment the sliding window.\n");
    fprintf(stderr, "                               Default 10.\n");
    fprintf(stderr, "   -k INT                  Dimension to project samples into. Must be less than number of samples.\n");
    fprintf(stderr, "                               Default 2. Must be less than number of samples in VCF.\n");
    fprintf(stderr, "   --threads INT           Number of threads to use in computation.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   --similarity            Compute similarity between sets of points instead of dissimilarity.\n");
    fprintf(stderr, "   --global                Compute only the global set of points.\n");
    fprintf(stderr, "                               Takes presedence over windowing parameters.\n");
    fprintf(stderr, "   --coords STR            A n-by-k csv file containing user defined coordinates to perform Procrustes analysis.\n");
    fprintf(stderr, "   --pthresh DOUBLE        Transform and print coordinates of all windows less than or equal to threshold.\n");
    fprintf(stderr, "                               Default 0.\n");
    fprintf(stderr, "   --perms INT             The number of permutations to execute.\n");
    fprintf(stderr, "                               Default 10000. Permutation test does not execute if pthresh is 0.\n");
    fprintf(stderr, "   --tthresh DOUBLE        Transform and print coordinates of all windows greater than or equal to threshold.\n");
    fprintf(stderr, "                               Default 0.95. Disabled if pthresh is not 0.\n");
    fprintf(stderr, "   --[^]regions REGIONS    Include/exclude records from VCF defined by REGIONS.\n");
    fprintf(stderr, "   --maf DOUBLE            Drops VCF records with a MAF less than threshold.\n");
    fprintf(stderr, "                               Default 0.\n");
    fprintf(stderr, "   --afMissing DOUBLE      Drops VCF records with fraction of missing genotypes greater than or equal to threshold.\n");
    fprintf(stderr, "                               Default 1.\n");
    fprintf(stderr, "   --[^]ibs REGIONS        Saves IBS calculations for overlapping windows included/excluded defined by REGIONS.\n");
    fprintf(stderr, "   --[^]asd REGIONS        Saves ASD calculations for overlapping windows included/excluded defined by REGIONS.\n");
    fprintf(stderr, "   --long                  Prints calculations in long format instead of matrix form.\n");
    fprintf(stderr, "   --json                  Prints window information in JSON format instead of TXT.\n");
    fprintf(stderr, "Types:\n");
    fprintf(stderr, "   STR                     A string.\n");
    fprintf(stderr, "   INT                     A non-negative integer.\n");
    fprintf(stderr, "   DOUBLE                  A real number between 0 and 1, inclusive.\n");
    fprintf(stderr, "   REGIONS                 REGION,REGIONS | REGION\n");
    fprintf(stderr, "   REGION                  STR | STR:INT | STR:-INT | STR:INT- | STR:INT-INT\n");
}

void print_window_info(Window* window, int n, int k) {
    printf("Window Number: %d\n", window -> winNum);
    printf("Chromosome: %s\n", ks_str(window -> chromosome));
    printf("Window Number on Chromosome: %d\n", window -> winNumOnChrom);
    printf("Start Position: %d\n", window -> startLocus);
    printf("End Position: %d\n", window -> endLocus);
    printf("Number of Loci: %d\n", window -> numLoci);
    printf("X = \n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++)
            printf("%5lf\t", window -> X[i][j]);
        printf("\n");
    }
    printf("\n\n");
}

static ko_longopt_t long_options[] = {
    {"help",            ko_no_argument,       300},
    {"version",         ko_no_argument,       301},
    {"similarity",      ko_no_argument,       302},
    {"global",          ko_no_argument,       303},
    {"threads",         ko_optional_argument, 304},
    {"pthresh",         ko_optional_argument, 305},
    {"perms",           ko_optional_argument, 306},
    {"tthresh",         ko_optional_argument, 307},
    {"regions",         ko_optional_argument, 308},
    {"^regions",        ko_optional_argument, 309},
    {"maf",             ko_optional_argument, 310},
    {"afMissing",       ko_optional_argument, 311},
    {"ibs",             ko_optional_argument, 312},
    {"^ibs",            ko_optional_argument, 313},
    {"asd",             ko_optional_argument, 314},
    {"^asd",            ko_optional_argument, 315},
    {"long",            ko_optional_argument, 316},
    {"json",            ko_optional_argument, 317},
    {0, 0, 0}
};

int main (int argc, char *argv[]) {

    const char *opt_str = "i:o:h:w:s:k:v";
    ketopt_t options = KETOPT_INIT;
    int c;

    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':':
                fprintf(stderr, "Error! Option is missing an argument! Exiting ...\n");
                return 1;
            case '?':
                fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]);
                return 1;
            case 300:
                print_help();
                return 0;
            case 'v':
            case 301:
                fprintf(stderr, "Version 1.0 May 2024.\n");
                return 0;
        }
	}
	options = KETOPT_INIT;

    LodestarConfiguration lodestar_config = {NULL, NULL, NULL, NULL, 1, 100, 10, 2, 1, false, false, NULL, NULL, 0, 10000, 0.95, NULL, 0, 1, NULL, NULL, NULL, NULL, false, false};

    print_configuration(stderr, lodestar_config);

    return 0;

}