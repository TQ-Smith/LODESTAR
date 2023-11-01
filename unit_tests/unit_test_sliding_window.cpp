//
// File: SlidingWindow.cpp
// Date: 31 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit test the sliding window algorithm for ASD calculations.
//

// Used for parsing VCF files.
#include "../src/VCFParser.hpp"

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

}