//
// File: unit_test_sliding_window.cpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit test the sliding window, ASD calculations, and MDS.
//

// Used to read in VCF files.
#include "../src/VCFParser.hpp"

// Import sliding window agorithm.
#include "../src/SlidingWindow.hpp"

#include <cassert>

#include <iostream>

using namespace std;

int main() {

    cout << endl;
    cout << "-----------------" << endl;
    cout << "First Test Case:" << endl;
    cout << "-----------------" << endl;

    // Open our VCF file.
    cout << endl;
    cout << "Opening sample1.vcf.gz ..." << endl;
    VCFParser* parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << endl;

    // Print the number of samples.
    int n = parser -> getNumberOfSamples();
    cout << "The VCF file has " << n << " samples." << endl;
    cout << endl;

    cout << "We create a sliding window with haplotype size of 1 SNPs, a window size of 1 haplotypes, and an offset of 1 haplotypes." << endl;
    cout << "We find the windows to be:" << endl;
    list<window*>* windows = window_genome(parser, 1, 1, 1, n, 2, false);
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++){
        cout << (*it) -> chromosome << ": " << (*it) -> start_position << " to " << (*it) -> end_position << " with ";
        if ((*it) -> chromosome == "Global") {
            cout << (*it) -> num_loci << " HAPs." << endl;
        } else {
            cout << (*it) -> num_loci << " SNPs." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;
    cout << endl;

    cout << "-----------------" << endl;
    cout << "Second Test Case:" << endl;
    cout << "-----------------" << endl;
    cout << endl;
    cout << "Opening sample1.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << endl;

    // Print the number of samples.
    parser -> getNumberOfSamples();
    cout << "The VCF file has " << n << " samples." << endl;
    cout << endl;

    cout << "We create a sliding window with haplotype size of 2 SNPs, a window size of 2 haplotypes, and an offset of 2 haplotypes." << endl;
    cout << "We find the windows to be:" << endl;
    windows = window_genome(parser, 2, 2, 2, n, 2, false);
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++){
        cout << (*it) -> chromosome << ": " << (*it) -> start_position << " to " << (*it) -> end_position << " with ";
        if ((*it) -> chromosome == "Global") {
            cout << (*it) -> num_loci << " HAPs." << endl;
        } else {
            cout << (*it) -> num_loci << " SNPs." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;
    cout << endl;

    cout << "-----------------" << endl;
    cout << "Third Test Case:" << endl;
    cout << "-----------------" << endl;
    cout << endl;
    cout << "Opening sample1.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << endl;

    // Print the number of samples.
    parser -> getNumberOfSamples();
    cout << "The VCF file has " << n << " samples." << endl;
    cout << endl;

    cout << "We create a sliding window with haplotype size of 3 SNPs, a window size of 1 haplotypes, and an offset of 1 haplotypes." << endl;
    cout << "We find the windows to be:" << endl;
    windows = window_genome(parser, 3, 1, 1, n, 2, false);
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++){
        cout << (*it) -> chromosome << ": " << (*it) -> start_position << " to " << (*it) -> end_position << " with ";
        if ((*it) -> chromosome == "Global") {
            cout << (*it) -> num_loci << " HAPs." << endl;
        } else {
            cout << (*it) -> num_loci << " SNPs." << endl;
        } 
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;
    cout << endl;

    cout << "-----------------" << endl;
    cout << "Fourth Test Case:" << endl;
    cout << "-----------------" << endl;
    cout << endl;
    cout << "Opening sample1.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << endl;

    // Print the number of samples.
    parser -> getNumberOfSamples();
    cout << "The VCF file has " << n << " samples." << endl;
    cout << endl;

    cout << "We create a sliding window with haplotype size of 1 SNPs, a window size of 3 haplotypes, and an offset of 2 haplotypes." << endl;
    cout << "We find the windows to be:" << endl;
    windows = window_genome(parser, 1, 3, 2, n, 2, false);
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++){
        cout << (*it) -> chromosome << ": " << (*it) -> start_position << " to " << (*it) -> end_position << " with ";
        if ((*it) -> chromosome == "Global") {
            cout << (*it) -> num_loci << " HAPs." << endl;
        } else {
            cout << (*it) -> num_loci << " SNPs." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;
    cout << endl;

    cout << "The code to check the ASD calculates is commented out in SlidingWindow.cpp." << endl;
    cout << endl;

}