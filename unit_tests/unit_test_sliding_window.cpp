//
// File: SlidingWindow.cpp
// Date: 31 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit test the sliding window algorithm and ASD calculations.
//              Does not test MDS.
//

// Used for parsing VCF files.
#include "../src/VCFParser.hpp"

// We are testing the sliding window.
#include "../src/SlidingWindow.hpp"

// Used for printing.
#include <iostream>
using namespace std;

int main() {

    // NOTE: Must uncomment print statements in SlidingWindow.cpp to print
    //  each window's dissimilarity matrix.

    // I *think* this covers all major cases.

    // Test single SNP, non overlapping windows.
    cout << endl;
    cout << "First Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample1.vcf.gz ..." << endl;
    VCFParser* parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << "Opened unit_tests/sample1.vcf.gz ..." << endl;
    cout << endl;

    int hap_size = 1;
    int win_size = 1;
    int step_size = 1;
    int n = 3;
    int k = 2;
    bool useFastMap = false;

    cout << "We window the genome with haplotype size of 1 SNP, window size of 1 haplotype, and step size of 1 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    // Window the genome.
    list<window*>* windows = window_genome(parser, hap_size, win_size, step_size, n, k, useFastMap);
    
    // Print windows.
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        cout << "Window on " << (*it) -> chromosome << " from " << (*it) -> start_position << " to " << (*it) -> end_position << " with " << (*it) -> num_loci;
        if ((*it) -> chromosome == "Global") {
            cout << " haplotypes." << endl;
        } else {
            cout << " loci." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;

    // Test mini-haplotype, non-overlapping windows
    cout << endl;
    cout << "Second Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample1.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << "Opened unit_tests/sample1.vcf.gz ..." << endl;
    cout << endl;

    hap_size = 2;
    win_size = 2;
    step_size = 2;
    n = 3;
    k = 2;
    useFastMap = false;

    cout << "We window the genome with haplotype size of 2 SNP, window size of 2 haplotype, and step size of 2 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    windows = window_genome(parser, hap_size, win_size, step_size, n, k, useFastMap);
    
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        cout << "Window on " << (*it) -> chromosome << " from " << (*it) -> start_position << " to " << (*it) -> end_position << " with " << (*it) -> num_loci;
        if ((*it) -> chromosome == "Global") {
            cout << " haplotypes." << endl;
        } else {
            cout << " loci." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;

    // Test the case when a single SNP is in a mini-haplotype at the end of a chromosome.
    cout << endl;
    cout << "Third Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample1.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << "Opened unit_tests/sample1.vcf.gz ..." << endl;
    cout << endl;

    hap_size = 3;
    win_size = 1;
    step_size = 1;
    n = 3;
    k = 2;
    useFastMap = false;

    cout << "We window the genome with haplotype size of 3 SNP, window size of 1 haplotype, and step size of 1 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    windows = window_genome(parser, hap_size, win_size, step_size, n, k, useFastMap);
    
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        cout << "Window on " << (*it) -> chromosome << " from " << (*it) -> start_position << " to " << (*it) -> end_position << " with " << (*it) -> num_loci;
        if ((*it) -> chromosome == "Global") {
            cout << " haplotypes." << endl;
        } else {
            cout << " loci." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;


    // Test SNP, overlapping windows.
    cout << endl;
    cout << "Fourth Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample1.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << "Opened unit_tests/sample1.vcf.gz ..." << endl;
    cout << endl;

    hap_size = 1;
    win_size = 3;
    step_size = 2;
    n = 3;
    k = 2;
    useFastMap = false;

    cout << "We window the genome with haplotype size of 1 SNP, window size of 3 haplotype, and step size of 2 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    windows = window_genome(parser, hap_size, win_size, step_size, n, k, useFastMap);
    
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        cout << "Window on " << (*it) -> chromosome << " from " << (*it) -> start_position << " to " << (*it) -> end_position << " with " << (*it) -> num_loci;
        if ((*it) -> chromosome == "Global") {
            cout << " haplotypes." << endl;
        } else {
            cout << " loci." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;

    // Test the case when different samples share identical genotypes
    //  for mini-haplotypes.
    cout << endl;
    cout << "Fifth Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample2.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample2.vcf.gz");
    cout << "Opened unit_tests/sample2.vcf.gz ..." << endl;
    cout << endl;

    hap_size = 5;
    win_size = 1;
    step_size = 1;
    n = 3;
    k = 2;
    useFastMap = false;

    cout << "We window the genome with haplotype size of 5 SNP, window size of 1 haplotype, and step size of 1 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    windows = window_genome(parser, hap_size, win_size, step_size, n, k, useFastMap);
    
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        cout << "Window on " << (*it) -> chromosome << " from " << (*it) -> start_position << " to " << (*it) -> end_position << " with " << (*it) -> num_loci;
        if ((*it) -> chromosome == "Global") {
            cout << " haplotypes." << endl;
        } else {
            cout << " loci." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;

    // Test the case of mini-haplotypes in overlapping windows.
    cout << endl;
    cout << "Sixth Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample3.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample3.vcf.gz");
    cout << "Opened unit_tests/sample3.vcf.gz ..." << endl;
    cout << endl;

    hap_size = 2;
    win_size = 2;
    step_size = 1;
    n = 3;
    k = 2;
    useFastMap = false;

    cout << "We window the genome with haplotype size of 2 SNP, window size of 2 haplotype, and step size of 1 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    windows = window_genome(parser, hap_size, win_size, step_size, n, k, useFastMap);
    
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        cout << "Window on " << (*it) -> chromosome << " from " << (*it) -> start_position << " to " << (*it) -> end_position << " with " << (*it) -> num_loci;
        if ((*it) -> chromosome == "Global") {
            cout << " haplotypes." << endl;
        } else {
            cout << " loci." << endl;
        }
        destroy_window(&*it, n);
    }
    delete parser;
    delete windows;

}