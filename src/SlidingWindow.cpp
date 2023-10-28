//
// File: SlidingWindow.cpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines the main algorithm for a sliding window along the genome.
//

#include "SlidingWindow.hpp"

// Use classical MDS and FastMap.
#include "../utils/LinearAlgebra/MultidimensionalScaling.hpp"

// Used to create and destroy matrices.
#include "../utils/LinearAlgebra/MatrixOperations.hpp"

// Used for the ceiling function.
#include <cmath>

// We define a macro to calculate the number of different alleles between two diploid individuals
//  at a single locus with the "Genotype" encoding from VCFParser.
#define LOCUS_ASD(i, j) (0.5 * (!!(i ^ j) + !(i & j)))

// We define a macro to convert ASD to ASD over the whole haplotype. If two samples have an 
//  ASD greater than or equal to 1 over a set of loci, then the haplotype ASD between the two is 1.
#define HAPLOTYPE_ASD(d) (d >= 1 ? 1 : d);

#include <iostream>
using namespace std;

window* create_mds_window(
    string chromosome, int start_position, int end_position, int num_loci, 
    double** D, int* maxDimReached, bool* doesConverge, double* d, double* e, 
    int n, int k, bool useFastMap) {
    
    double** points = create_real_matrix(n, k);
    if (useFastMap) {
        compute_fastmap(D, points, maxDimReached, n, k);
    } else {
        compute_classical_mds(D, points, d, e, doesConverge, n, k);
    }
    window* w = new window;
    w -> chromosome = chromosome;
    w -> start_position = start_position;
    w -> end_position = end_position;
    w -> num_loci = num_loci;
    w -> points = points;
    return w;

}

list<window*>* window_genome(VCFParser* parser, int hap_size, int window_hap_size, int offset_hap_size, int n, int k, bool useFastMap) {
    
    list<window*>* windows = new list<window*>;

    Genotype* genotypes = new Genotype[n];

    double** window_asd = new double*[n];
    double** auxilary_asd = new double*[n];
    double** next_asd = new double*[n];
    double** global_asd = new double*[n];
    for (int i = 0; i < n; i++) {
        window_asd[i] = new double[n]; auxilary_asd[i] = new double[n]; next_asd[i] = new double[n]; global_asd[i] = new double[n];
        window_asd[i][i] = auxilary_asd[i][i] = next_asd[i][i] = global_asd[i][i] = 0;
    }

    string chrom;
    int pos;
    bool isMonomorphic, isComplete, nextRecord;

    int q;
    cout << ((long) &q) % 16 << endl;

    double* d = new double[n];
    double* e = new double[n];
    bool doesConverge;
    int maxDimReached;

    int p;
    cout << ((long) &p) % 16 << endl;

    string prev_chrom = "";
    int start_pos, next_start_pos, prev_pos, num_haps, num_loci = 0, num_global_haps = 0;
    double hap_asd;

    while ( true ) {

        nextRecord = parser -> getNextLocus(&chrom, &pos, &isMonomorphic, &isComplete, genotypes);

        if (isMonomorphic || !isComplete) {
            continue;
        }

        if (num_loci != 0 && (!nextRecord || chrom != prev_chrom || num_loci == (hap_size * window_hap_size))) {

            num_haps = ceil(num_loci / hap_size);
            if (!nextRecord) {
                num_global_haps += num_haps;
            } else {
                num_global_haps += offset_hap_size + ceil((num_loci - hap_size * offset_hap_size) / hap_size);
            }

            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    hap_asd = HAPLOTYPE_ASD(auxilary_asd[i][j]);
                    next_asd[i][j] = next_asd[j][i] = window_asd[i][j] - auxilary_asd[j][i];
                    window_asd[i][j] = window_asd[j][i] = (window_asd[i][j] + hap_asd) / num_haps;
                    if (!nextRecord) {
                        global_asd[i][j] = global_asd[j][i] = (global_asd[i][j] + hap_asd) / num_global_haps;
                    }
                    auxilary_asd[j][i] = 0;
                }
            }

            windows -> push_back(create_mds_window(prev_chrom, start_pos, prev_pos, num_loci, window_asd, &maxDimReached, &doesConverge, d, e, n, k, useFastMap));

            double** temp = window_asd;
            window_asd = next_asd;
            next_asd = temp;

            num_loci -= hap_size * offset_hap_size;
            start_pos = (window_hap_size == offset_hap_size) ? pos : next_start_pos;

        }

        if (!nextRecord) {
            break;
        }

        if (prev_chrom != chrom) {
            num_loci = 0;
            next_start_pos = pos;
            start_pos = pos;
        }

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (prev_chrom != chrom) {
                    window_asd[i][j] = window_asd[j][i] = auxilary_asd[i][j] = auxilary_asd[j][i] = next_asd[i][j] = next_asd[j][i] = 0;
                }
                if (num_loci != 0 && num_loci % hap_size == 0) {
                    hap_asd = HAPLOTYPE_ASD(auxilary_asd[i][j]);
                    window_asd[i][j] = window_asd[j][i] += hap_asd;
                    global_asd[i][j] = global_asd[j][i] += hap_asd;
                    if (num_loci <= hap_size * offset_hap_size) {
                        auxilary_asd[j][i] += hap_asd;
                    }
                    auxilary_asd[i][j] = 0;
                }
                auxilary_asd[i][j] += LOCUS_ASD(genotypes[i], genotypes[j]);
            }
        }

        prev_chrom = chrom;
        prev_pos = pos;
        num_loci++;

    }

    windows -> push_back(create_mds_window("Global", windows -> front() -> start_position, prev_pos, num_global_haps, global_asd, &maxDimReached, &doesConverge, d, e, n, k, useFastMap));

    delete [] genotypes;
    destroy_real_matrix(window_asd, n);
    destroy_real_matrix(auxilary_asd, n);
    destroy_real_matrix(next_asd, n);
    destroy_real_matrix(global_asd, n);
    delete [] d;
    delete [] e;

    return windows;

}

void destroy_window(window** w, int n) {

    // Deallocate points matrix.
    destroy_real_matrix((*w) -> points, n);

    // Deallocate structure.
    delete *w;

    // Set pointer to NULL.
    *w = NULL;

}