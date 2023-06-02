//
// File: Engine.cpp
// Started: 17 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains the main logic/flow of the analysis.
//          Also handles threading.
//

// Used to define main function operations.
#include "Engine.hpp"

// Used for allocating new matrices and deep-copies.
#include "MatrixOperations.hpp"

// Used to preform the MDS operation.
#include "MultidimensionalScaling.hpp"

// Used to preform procurstes analysis and permutation tests.
#include "Procrustes.hpp"

// Used to parse genotype data file.
#include "GenotypeFileParser.hpp"

// Used to hold pointers to windows.
#include <list>

// Used for printing.
#include <iostream>

// A private method to convert a similarity count matrix to dissimilarity.
// Accepts:
//  Window* window -> The window holds the similarity count matrix and number of loci in the window.
//  int n -> The dimension of the matrix.
// Returns: void.
void convert_to_dissimilarity(Window* window, int n) {
    for (int i = 0; i < n; i++) {
        window -> points[i][i] = 1;
        for (int j = i + 1; j < n; j++) {
            // Note, we are only working with diploid individuals.
            window -> points[i][j] = window -> points[j][i] = 1 - (window -> points[i][j]) / (2 * (window -> num_loci));
        }
    }
}

// A private method to calculate the allele similarity count between two genotypes.
//  Pass by reference is not needed, but I'll take the speed up.
// Accepts:
//  Genotype& sample1 -> The genotype of the first sample.
//  Genotype& sample2 -> The genotype of the second sample.
// Returns: int, The number of similar alleles. Either 0, 1, or 2.
int calculate_ads(Genotype& sample1, Genotype& sample2) {
    if (sample1.chr1 == sample2.chr1) {
        if (sample1.chr2 == sample2.chr2) {
            return 2;
        }
        return 1;
    }
    if (sample1.chr1 == sample2.chr2) {
        if (sample1.chr2 == sample2.chr1) {
            return 2;
        }
        return 1;
    }
    if (sample1.chr2 == sample2.chr2) {
        if (sample1.chr1 == sample2.chr2) {
            return 2;
        }
        return 1;
    }
    return 0;
}

// A private method to deallocate a window's memory.
// Accepts:
//  Window* window -> The pointer to the window structure to destroy.
//  int n -> The number of rows.
//  int k -> The number of columns.
// Returns: void.
void destroy_window(Window* window, int n, int k) {
    destroy_real_matrix(window -> points, n, k);
    delete window;
}

void print_window(Window* window, int n, int k) {
    cout << "Chromosome: " << window -> chromosome << endl;
    cout << "Start Position: " << window -> start_position << endl;
    cout << "End Position: " << window -> end_position << endl;
    cout << "Number of Loci: " << window -> num_loci << endl;
    cout << "Dissimilarity Matrix:" << endl;
    print_real_matrix(window -> points, n, k);  
}

void window_genome(ifstream& in_file, list<Window*>& windows, string unit, int window_width, int window_offset, int n, int k) {
    
    // Create global window.
    //  Initialize to empty count matrix and no loci.
    Window* global = new Window;
    global -> chromosome = "Global";
    global -> points = create_and_fill_real_matrix(0, n, n);
    global -> num_loci = 0;

    // Variables used to read in a locus
    //  The chromosome of the locus.
    //  The position of the locus.
    //  Our flag variable to indicate end of file.
    //  Our flag variable to indiciate incomplete genotypes.
    //  Our flag varaible to indicate that all genotypes were not homogenous.
    //  Our array of geneotypes for each individual.
    string chromosome = "";
    int position;
    bool isComplete;
    bool isEOF;
    Genotype* genotypes = new Genotype[n];

    // We create a window to keep track of the allele similarity counts
    //  within the current window width.
    Window* width_window = new Window;
    width_window -> points = create_real_matrix(n, n);
    width_window -> num_loci = 0;

    // We have a window to keep track of the allele similarity counts
    //  that is within the current window frame and the width of the 
    //  offset. This way, in the sliding window, we can save computation
    //  by subtracting the overlap from the current window counts to 
    //  contribute to the next window, as the window advances.
    Window* overlap_window = new Window;
    overlap_window -> points = create_real_matrix(n, n);
    overlap_window -> num_loci = 0;

    // We keep track of the previous chromosome and position. If a new chromosome is encountered,
    //  then we have to create a new window.
    string previous_chromosome = "";
    int previous_position;

    // A temporary varaible used to calculate the distances between individuals.
    int similarity;

    // Used for MDS.
    double* d = new double[n];
    double* e = new double[n];

    int next_start = 0;

    // Our loop to read in each locus.
    while(true) {
        
        // First, we read in the locus.
        get_next_loci(in_file, chromosome, position, genotypes, isComplete, isEOF, n);

        // If the EOF flag was set, break out of the loop.
        if (isEOF) {
            break;
        }

        // If there is an incomplete genotype for an individual, then continue to the next locus.
        if (!isComplete) {
            continue;
        }

        if (previous_chromosome != chromosome || width_window -> num_loci == window_width) {

            if (width_window -> num_loci != 0) {
                width_window -> end_position = previous_position;
                width_window -> chromosome = previous_chromosome;
                subtract_matrices(overlap_window -> points, width_window -> points, overlap_window -> points, n, n);
                overlap_window -> num_loci = width_window -> num_loci - overlap_window -> num_loci;
                convert_to_dissimilarity(width_window, n);
                // print_real_matrix(width_window -> points, n, n);
                compute_classical_mds(width_window -> points, d, e, n, k);
                windows.push_back(width_window);
                width_window = overlap_window;
                width_window -> start_position = next_start;
                overlap_window = new Window;
                overlap_window -> num_loci = 0;
                overlap_window -> start_position = position;
                overlap_window -> points = create_and_fill_real_matrix(0, n, n);
            }

            if (previous_chromosome != chromosome) {
                width_window -> start_position = position;
                width_window -> num_loci = 0;
                overlap_window -> start_position = position;
                overlap_window -> num_loci = 0;
                for (int i = 0; i < n; i++) {
                    for (int j = i; j < n; j++) {
                        width_window -> points[i][j] = width_window -> points[j][i] = overlap_window -> points[i][j] = overlap_window -> points[j][i] = 0;
                    }
                }
            }

        }
        
        // If the locus was in the overlap, then we add the distances to the count matrix.
        //  Otherwise, we do not add it to the overlap count.
        if (overlap_window -> num_loci < window_offset) {
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    // Calculate allele similarity.
                    similarity = calculate_ads(genotypes[i], genotypes[j]);
                    overlap_window -> points[i][j] += similarity;
                    width_window -> points[i][j] += similarity;
                    global -> points[i][j] += similarity;
                    overlap_window -> points[j][i] += similarity;
                    width_window -> points[j][i] += similarity;
                    global -> points[j][i] += similarity;
                }
            }
            overlap_window -> num_loci++;
            overlap_window -> end_position = position;
        } else {
            if (width_window -> num_loci == (window_width - window_offset) - 1 ) {
                next_start = position;
            }
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    // Calculate allele similarity.
                    similarity = calculate_ads(genotypes[i], genotypes[j]);
                    width_window -> points[i][j] += similarity;
                    global -> points[i][j] += similarity;
                    width_window -> points[j][i] += similarity;
                    global -> points[j][i] += similarity;
                }
            }
        }

        previous_chromosome = chromosome;
        previous_position = position;
        width_window -> num_loci++;
        global -> num_loci++;
    }

    // If the end of the file was reached with an unfinished window.
    //  Do the same as above, but no need to advance the slider.
    if (width_window -> num_loci != 0) {
        width_window -> chromosome = chromosome;
        width_window -> end_position = position;
        convert_to_dissimilarity(width_window, n);
        // print_real_matrix(width_window -> points, n, n);
        compute_classical_mds(width_window -> points, d, e, n, k);
        windows.push_back(width_window);
    }

    // Global counts matrix is converted to dissimilarity 
    //  and pushed onto the end of the list.
    //  We set the first and last starting position here.
    global -> start_position = windows.front() -> start_position;
    global -> end_position = windows.back() -> end_position;
    convert_to_dissimilarity(global, n);
    // print_real_matrix(global -> points, n, n);
    compute_classical_mds(global -> points, d, e, n, k);
    windows.push_back(global);

    // Destroy overlap window and geneotypes array.
    destroy_window(overlap_window, n, k);
    delete[] genotypes;
    delete[] d;
    delete[] e;

}

void lodestar_pipeline(string input_file_name, string unit, int window_width, int window_offset, int k) {

    cout << endl;
    cout << "Opening input file " << input_file_name << endl;
    ifstream in_file;
    open_file(in_file, input_file_name);
    cout << "Opened input file " << input_file_name << endl;
    cout << endl;

    string* names;
    int n;

    cout << "Getting sample names." << endl;
    get_sample_names(in_file, names, n);
    cout << "There are " << n << " samples in the file." << endl;
    cout << "The samples are named:" << endl;
    for (int i = 0; i < n; i++) {
        cout << names[i] << endl;
    }
    cout << endl;
    
    list<Window*> windows;

    cout << "Now, we window the genome." << endl;
    window_genome(in_file, windows, unit, window_width, window_offset, n, k);
    cout << "We finished windowing the genome." << endl;

    cout << endl;
    Window* temp;
    int size = windows.size();
    for (int i = 0; i < size; i++) {
        temp = windows.front();
        print_window(temp, n, k);
        destroy_window(temp, n, n);
        windows.pop_front();
        cout << endl;
    }

    cout << "Closing input file." << endl;
    close_file(in_file);
    cout << "Closed input file." << endl;
    cout << endl;
    

}