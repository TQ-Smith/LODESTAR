//
// File: unit_test_sliding_window.cpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit test all of the cases that could be encountered
//              for a sliding window. NOTE: Does not test MDS and 
//              Procrustes calculations. 
//

// Used to read in VCF files.
#include "../src/VCFParser.hpp"

// Import sliding window agorithm.
#include "../src/SlidingWindow.hpp"

#include <cassert>

#include <iostream>

using namespace std;

int main() {

    // Open our VCF file.
    cout << endl;
    cout << "Opening sample2.vcf ..." << endl;
    VCFParser* parser = new VCFParser("unit_tests/sample2.vcf");
    cout << endl;

    // Print the number of samples.
    int n = parser -> getNumberOfSamples();
    cout << "The VCF file has " << n << " samples." << endl;
    cout << endl;

    cout << "We create a sliding window with haplotype size of 1 SNPs, a window size of 2 haplotypes, and an offset of 1 haplotypes." << endl;
    cout << "We find the windows to be:" << endl;
    list<window*>* windows = sliding_window(parser, 1, 2, 1, n);
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++){
        cout << (*it) -> chromosome << ": " << (*it) -> start_position << " to " << (*it) -> end_position << " with " << (*it) -> num_loci << " SNPs." << endl;
        destroy_window(&*it);
    }

}