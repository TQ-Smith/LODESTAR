//
// File: SlidingWindow.cpp
// Date: 31 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit test the sliding window algorithm for ASD calculations.
//

// Used for parsing VCF files.
#include "../src/VCFParser.hpp"

// We are testing the sliding window.
#include "../src/SlidingWindow.hpp"

// Used for printing.
#include <iostream>
using namespace std;

int main() {

    cout << endl;
    cout << "First Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample1.vcf.gz ..." << endl;
    VCFParser* parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << "Opened unit_tests/sample1.vcf.gz ..." << endl;
    cout << endl;

    int hap_size = 1;
    int window_hap_size = 1;
    int offset_hap_size = 1;
    int n = 3;
    int k = 3;
    bool useFastMap = false;

    cout << "We window the genome with haplotype size of 1 SNP, window size of 1 haplotype, and offset size of 1 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    list<window*>* windows = window_genome(parser, hap_size, window_hap_size, offset_hap_size, n, k, useFastMap);
    
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

    cout << endl;
    cout << endl;
    cout << "Second Test:" << endl;
    cout << "-----------" << endl;
    cout << endl;
    cout << "Opening unit_tests/sample1.vcf.gz ..." << endl;
    parser = new VCFParser("unit_tests/sample1.vcf.gz");
    cout << "Opened unit_tests/sample1.vcf.gz ..." << endl;
    cout << endl;

    hap_size = 2;
    window_hap_size = 2;
    offset_hap_size = 2;
    n = 3;
    k = 3;
    useFastMap = false;

    cout << "We window the genome with haplotype size of 2 SNP, window size of 2 haplotype, and offset size of 2 haplotype." << endl;
    cout << "The generated windows:" << endl << endl;

    windows = window_genome(parser, hap_size, window_hap_size, offset_hap_size, n, k, useFastMap);
    
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