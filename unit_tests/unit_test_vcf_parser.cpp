//
// File: unit_test_vcf_parser.cpp
// Date: 5 September 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit test the VCF file parser.
//

// Import genotype operations.
#include "../src/Genotype.hpp"

// Import the parser.
#include "../src/VCFParser.hpp"

// Used for printing and file name.
#include <iostream>

using namespace std;

// Used for testing.
#include <cassert>

int main() {

    // Test if is open works.
    VCFParser foo("foo");
    assert(!foo.isOpen());

    // Load sample file for testing.
    //  We are also testing gzstream.
    VCFParser parser("unit_tests/sample1.vcf.gz");
    assert(parser.isOpen());

    // First, we test the constructor and the getter for
    //  the number of samples and the names of the samples.
    int n = parser.getNumberOfSamples();
    assert(n == 3);

    string* samples = parser.getSampleNames();
    assert(samples[0] == "NA00001");
    assert(samples[1] == "NA00002");
    assert(samples[2] == "NA00003");
    delete[] samples;

    // Create all other necessary variable for reading in the file.
    string chromosome;
    int position;
    bool isMonomorphic;
    bool isComplete;

    // Create our genotype array.
    Genotype* genotypes = new Genotype[n];

    // Read in the first locus.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(chromosome == "1");
    assert(position == 100);
    assert(!isMonomorphic);
    assert(isComplete);
    assert(genotypes[0] == 3);
    // Since the genotype is 1/0, the most significant bit will be set.
    assert(genotypes[1] == 0x80000003);
    // Since the genotype is 1/1, the most significant bit will be set.
    assert(genotypes[2] == 0x80000002);

    // Read in the second locus.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(!isMonomorphic);
    assert(isComplete);

    // Read in the third locus.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(!isMonomorphic);
    assert(isComplete);

    // Read in the fourth locus.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(!isMonomorphic);
    assert(isComplete);

    // Read in the fifth locus.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(!isMonomorphic);
    assert(isComplete);

    // Read in the sixth locus. This one is monomorphic.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(isMonomorphic);
    assert(isComplete);

    // Read in the seventh locus.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(!isMonomorphic);
    assert(isComplete);

    // Read in the eighth locus.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(isComplete);

    // Read in the ninth locus. This one is not complete.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(!isComplete);

    // Read in the tenth locus. This one is not complete.
    assert(parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));
    assert(!isComplete);

    // End of file.
    assert(!parser.getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes));

    delete[] genotypes;

}