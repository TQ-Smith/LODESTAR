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

void print_window(Window* window, int n) {
    cout << "Chromosome: " << window -> chromosome << endl;
    cout << "Start Position: " << window -> start_position << endl;
    cout << "End Position: " << window -> end_position << endl;
    cout << "Number of Loci: " << window -> num_loci << endl;
    cout << "Dissimilarity Matrix:" << endl;
    print_real_matrix(window -> points, n, n);  
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

    // A counter for the number of windows. Used to keep track of window boundds for bp.
    int num_windows = -1;

    // The unit changes two things. How we determine if a locus is in the overlap
    //  and if a locus is in the window's width. We create a variable for each test.
    bool isInOverlap;
    bool isInWindow;

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

        // Determine if the locus was in in the window.
        //  We account for the chromosomes not being equal too.
        //  Default is always bp.
        //  Add in centimorgans later.
        if (unit == "snp") {
            isInWindow = previous_chromosome == chromosome && ((width_window -> num_loci) < window_width);
        } else {
            isInWindow = previous_chromosome == chromosome && (position < (window_offset * num_windows + window_width));
        }

        // Next, we determine if the loci starts in a new window.
        if (!isInWindow) {
            // Consider the event of an empty window.
            //  If the window is empty, then so will the offset,
            //  since the offset is smaller than the window width.
            if ((width_window -> num_loci) == 0) {
                
                // We zero out the counts of the overlap and width counts.
                for (int i = 0; i < n; i++) {
                    for (int j = i; j < n; j++) {
                        width_window -> points[i][j] = width_window -> points[j][i] = 0;
                        overlap_window -> points[i][j] = overlap_window -> points[j][i] = 0;
                    }
                }

                // Set the chromosome and starting position.
                width_window -> start_position = position;
                width_window -> num_loci = 0;
            
            // If the window is not empty, then we must process it and start a new window.
            } else {

                // Set the ending position.
                //  Corresponds to the last locus in the window.
                width_window -> end_position = previous_position;

                // NOTE: If the windows are non-overlapping, then we could
                //  have just created a new matrix, but just like subtraction
                //  that is a quadratic operation. Thus, we do not consider
                //  this case.

                // We now subtract the overlap from the width_window.
                //  This prevents us from preforming double computations
                //  and having to save loci in memory.
                subtract_matrices(overlap_window -> points, width_window -> points, overlap_window -> points, n, n);

                 // Adjust number of loci in the window but not the overlap.
                overlap_window -> num_loci = width_window -> num_loci - overlap_window -> num_loci;
                
                // Convert the width_window to dissimilarity and save to list.
                convert_to_dissimilarity(width_window, n);

                windows.push_back(width_window);

                // We advance the sliding window by setting the new width_window to
                //  the old overlap window. We then allocate a new window for the
                //  width_window.
                width_window = overlap_window;
                overlap_window = new Window;

                // Reset the start position to the new window.
                width_window -> start_position = position;

                // There are no loci in the overlap window yet.
                overlap_window -> num_loci = 0;

                // Lastly, we create the new matrix for the new window overlap window.
                overlap_window -> points = create_real_matrix(n ,n);

            }

            // We created a new window.
            num_windows++;

            // Declare chromsome of window.
            width_window -> chromosome = chromosome;
        }

        // The window is setup. We process the locus.

        // Determine if locus is in the overlap.
        //  Default bp.
        //  Add centimorgans later.
        if (unit == "snp") {
            isInOverlap = (overlap_window -> num_loci) < window_offset;
        } else {
            isInOverlap = position < (window_offset * (num_windows + 1));
        }

        // If the locus was in the overlap, then we add the distances to the count matrix.
        //  Otherwise, we do not add it to the overlap count.
        if (isInOverlap) {
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    // Calculate allele similarity.
                    similarity = calculate_ads(genotypes[i], genotypes[j]);
                    cout << genotypes[i].chr1 << "," << genotypes[i].chr2 << " " << genotypes[j].chr1 << "," << genotypes[j].chr2 << ": " << similarity << endl;
                    overlap_window -> points[i][j] = width_window -> points[i][j] = global -> points[i][j] = similarity;
                    overlap_window -> points[j][i] = width_window -> points[j][i] = global -> points[j][i] = similarity;
                }
            }
            overlap_window -> num_loci++;
        } else {
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    // Calculate allele similarity.
                    similarity = calculate_ads(genotypes[i], genotypes[j]);
                    width_window -> points[i][j] = global -> points[i][j] = similarity;
                    width_window -> points[j][i] = global -> points[j][i] = similarity;
                }
            }
        }

        // We processed the window.
        //  Increment the number of loci and
        //  update the previous chromosome.
        width_window -> num_loci++;
        global -> num_loci++;
        previous_chromosome = chromosome;
        previous_position = position;
    }

    // If the end of the file was reached with an unfinished window.
    //  Do the same as above, but no need to advance the slider.
    if (width_window -> num_loci != 0) {
        width_window -> end_position = position;
        convert_to_dissimilarity(width_window, n);
        windows.push_back(width_window);
    }

    // Global counts matrix is converted to dissimilarity 
    //  and pushed onto the end of the list.
    //  We set the first and last starting position here.
    global -> start_position = windows.front() -> start_position;
    global -> end_position = windows.back() -> end_position;
    convert_to_dissimilarity(global, n);
    windows.push_back(global);

    // Destroy overlap window and geneotypes array.
    destroy_window(overlap_window, n, k);
    delete[] genotypes;

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
        print_window(temp, n);
        destroy_window(temp, n, k);
        windows.pop_front();
        cout << endl;
    }

    cout << "Closing input file." << endl;
    close_file(in_file);
    cout << "Closed input file." << endl;
    cout << endl;
    

}