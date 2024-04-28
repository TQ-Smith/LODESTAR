
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include "ProcrustesAnalysis.h"

#include "../lib/ketopt.h"


void print_help() {
    fprintf(stderr, "LODESTAR v1.0 May 2024\n");
    fprintf(stderr, "----------------------\n\n");
    fprintf(stderr, "Written by T. Quinn Smith\n");
    fprintf(stderr, "Principal Investigator: Zachary A. Szpiech\n");
    fprintf(stderr, "The Pennsylvania State University\n\n");
    fprintf(stderr, "Usage: lodestar [options] -i <input.vcf.gz> -o <output_basename>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "   --help              Prints help menu and exits.\n");
    fprintf(stderr, "   --version           Prints version number and exits.\n");
    fprintf(stderr, "   -i STR              Path to input vcf file.\n");
    fprintf(stderr, "   -o STR              Output base names for files.\n");
    fprintf(stderr, "   -h INT              Number of loci in a haplotype. Default 10.\n");
    fprintf(stderr, "   -w INT              Number of haplotypes in a window. Default 100.\n");
    fprintf(stderr, "   -s INT              Number of haplotypes to increment the sliding window. Default 10.\n");
    fprintf(stderr, "   -k INT              Dimension to project samples into. Must be less than number of samples. Default 2.\n");
    fprintf(stderr, "   --threads INT       Number of threads to use in computation.\n");
    fprintf(stderr, "   --similarity        Compute similarity between sets of points instead of dissimilarity.\n");
    fprintf(stderr, "   --global            Compute only the global set of points.\n");

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
    {"perms",           ko_optional_argument, 305},
    {0, 0, 0}
};

int main (int argc, char *argv[]) {

    bool help = false;
    char* input_file = NULL;
    char* output_basename = NULL;
    int HAP_SIZE = 10;
    int WINDOW_SIZE = 100;
    int STEP_SIZE = 10;
    int k = 2;
    bool similarity = false;
    bool global = false;
    int threads = 1;

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

    /*
    int k = 1;
    int NUM_THREADS = 1;
    int HAP_SIZE = 1, STEP_SIZE = 3, WINDOW_SIZE = 4;

    VCFLocusParser* parser = init_vcf_locus_parser("./data/sliding_window_test2.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> numSamples);

    int numWindows;
    Window** windows = sliding_window(parser, encoder, k, HAP_SIZE, STEP_SIZE, WINDOW_SIZE, NUM_THREADS, &numWindows);

    for (int i = 0; i < numWindows; i++)
        print_window_info(windows[i], encoder -> numSamples, k);

    for (int i = 0; i < numWindows; i++)
        destroy_window(windows[i], parser -> numSamples);
    free(windows);
    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);
    */

    return 0;

}