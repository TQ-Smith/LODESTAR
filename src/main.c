
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include "ProcrustesAnalysis.h"

#include "../lib/ketopt.h"

typedef struct {
    char* input_filename;
    VCFLocusParser* parser;
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
    char* coords_filename;
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
    fprintf(output, "P-value threshold: %lf\n", lodestar_config.pthresh);
    fprintf(output, "Number of permutations: %d\n", lodestar_config.num_perms);
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
    if (lodestar_config.coords_filename != NULL)
        fprintf(output, "File of coordinates to perform Procrustes analysis against: %s\n", lodestar_config.coords_filename);
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
    if (lodestar_config.coords_file != NULL)
        fclose(lodestar_config.coords_file);
    if (lodestar_config.ibs_regions != NULL)
        destroy_region_filter(lodestar_config.ibs_regions);
    if (lodestar_config.asd_regions != NULL)
        destroy_region_filter(lodestar_config.asd_regions);
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
    {"[^]regions",      ko_optional_argument, 308},
    {"maf",             ko_optional_argument, 309},
    {"afMissing",       ko_optional_argument, 310},
    {"[^]ibs",          ko_optional_argument, 311},
    {"[^]asd",          ko_optional_argument, 312},
    {"long",            ko_optional_argument, 313},
    {"json",            ko_optional_argument, 314},
    {"coords",          ko_optional_argument, 315},
    {0, 0, 0}
};

int main (int argc, char *argv[]) {

    const char *opt_str = "i:o:h:w:s:k:v";
    ketopt_t options = KETOPT_INIT;
    int c;

    while ((c = ketopt(&options, argc, argv, 1, opt_str, long_options)) >= 0) {
        switch (c) {
            case ':': fprintf(stderr, "Error! Option is missing an argument! Exiting ...\n"); return 1;
            case '?': fprintf(stderr, "Error! \"%s\" is unknown! Exiting ...\n", argv[options.i - 1]); return 1;
            case 300: print_help(); return 0;
            case 'v': case 301: fprintf(stderr, "Version 1.0 May 2024.\n"); return 0;
        }
	}
	options = KETOPT_INIT;

    LodestarConfiguration lodestar_config = {NULL, NULL, NULL, NULL, NULL, 1, 100, 10, 2, 1, false, false, NULL, NULL, 0, 10000, 0.95, NULL, 0, 1, NULL, NULL, NULL, NULL, false, false};

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
            case 315: lodestar_config.coords_filename = options.arg; break;
        }
	}

    print_configuration(stderr, lodestar_config);

    return 0;

}