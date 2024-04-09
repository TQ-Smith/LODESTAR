
#include "VCFLocusParser.h"

#include "HaplotypeEncoder.h"

#include "SlidingWindow.h"

#include "MultidimensionalScaling.h"

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
    int NUM_THREADS = 2;
    int HAP_SIZE = 1, STEP_SIZE = 1, WINDOW_SIZE = 3;

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

    int N = 4;
    int k = 2;
    int v;
    RealSymEigen* eigen = init_real_sym_eigen(N);
    double* D = malloc(PACKED_SIZE(N) * sizeof(double));
    printf("D:\n");
    for (int i = 0; i < N; i++) {
        v = 1;
        for (int j = i; j < N; j++) {
            D[INDEX(i, j)] = v++;
            printf("%lf\t", D[INDEX(i, j)]);
        }
        printf("\n");
    }

    double** X = create_matrix(double, N, k);

    compute_classical_mds(eigen, D, k, X);

    printf("X:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            printf("%3.2f\t", X[i][j]);
        }
        printf("\n");
    }

    destroy_matrix(double, X, N);
    destroy_real_sym_eigen(eigen);

    return 0;

}