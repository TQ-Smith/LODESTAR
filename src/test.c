
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include "MultidimensionalScaling.h"

#include "ProcrustesAnalysis.h"

#include <stdio.h>

#include <time.h>


int main() {

    int k = 2;
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
    

    /*
    int n = 3;
    int k = 2;
    bool transform = false;
    bool similarity = false;
    RealSymEigen* eigen = init_real_sym_eigen(k);
    double** X = create_matrix(double, n, k);
    double* x0 = (double*) calloc(k, sizeof(double));
    X[0][0] = 1; X[0][1] = 2; X[1][0] = 3; X[1][1] = 4; X[2][0] = 5; X[2][1] = 6;
    x0[0] = 1; x0[1] = 2;
    double** Y = create_matrix(double, n, k);
    double* y0 = (double*) calloc(k, sizeof(double));
    Y[0][0] = 3; Y[0][1] = 5; Y[1][0] = 7; Y[1][1] = 11; Y[2][0] = 13; Y[2][1] = 15;
    y0[0] = 5; y0[1] = 10;

    double t = procrustes_statistic(X, x0, Y, y0, eigen, n, k, transform, similarity);

    
    printf("Transformed X:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            printf("%lf\t", X[i][j]);
        }
        printf("\n");
    }
    printf("\nDissimilarity: %lf\n\n", t);
    

    double** shuffleX = create_matrix(double, n, k);

    srand(time(NULL));

    double p = permutation_test(X, Y, shuffleX, eigen, n, k, similarity, t, 10000);

    printf("p-value: %lf\n", p);

    destroy_real_sym_eigen(eigen);
    destroy_matrix(double, X, n);
    destroy_matrix(double, Y, n);
    destroy_matrix(double, shuffleX, n);
    free(x0);
    free(y0);
    */

    return 0;

}