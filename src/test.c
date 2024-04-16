
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include "MultidimensionalScaling.h"

#include "ProcrustesAnalysis.h"

#include <stdio.h>

#include "Matrix.h"
MATRIX_INIT(double, double)

void print_window_info(Window* window) {
    printf("Window Number: %d\n", window -> winNum);
    printf("Chromosome: %s\n", ks_str(window -> chromosome));
    printf("Window Number on Chromosome: %d\n", window -> winNumOnChrom);
    printf("Start Position: %d\n", window -> startLocus);
    printf("End Position: %d\n", window -> endLocus);
    printf("Number of Loci: %d\n", window -> numLoci);
    printf("\n");
}

int main() {

    /*
    int NUM_THREADS = 1;
    int HAP_SIZE = 1, STEP_SIZE = 1, WINDOW_SIZE = 2;

    VCFLocusParser* parser = init_vcf_locus_parser("./data/sliding_window_test.vcf.gz");
    HaplotypeEncoder* encoder = init_haplotype_encoder(parser -> numSamples);

    int numWindows;
    Window** windows = sliding_window(parser, encoder, HAP_SIZE, STEP_SIZE, WINDOW_SIZE, NUM_THREADS, &numWindows);

    for (int i = 0; i < numWindows; i++)
        print_window_info(windows[i]);

    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);
    for (int i = 0; i < numWindows; i++)
        destroy_window(windows[i]);
    free(windows);
    */

    int n = 3;
    int k = 2;
    bool transform = false;
    bool similarity = false;
    RealSymEigen* eigen = init_real_sym_eigen(k);
    double** X = create_matrix(double, n, k);
    X[0][0] = 1; X[0][1] = 1; X[1][0] = 2; X[1][1] = 3; X[2][0] = 4; X[2][1] = 7;
    double** Y = create_matrix(double, n, k);
    Y[0][0] = 2; Y[0][1] = 7; Y[1][0] = 1; Y[1][1] = 4; Y[2][0] = 1; Y[2][1] = 1;

    procrustes_statistic(X, NULL, Y, NULL, eigen, n, k, transform, similarity);

    destroy_real_sym_eigen(eigen);
    destroy_matrix(double, X, n);
    destroy_matrix(double, Y, n);

    return 0;

}