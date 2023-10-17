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

    cout << "Opening sample2.vcf ..." << endl;
    VCFParser* parser = new VCFParser("unit_tests/sample2.vcf");
    cout << endl;

}