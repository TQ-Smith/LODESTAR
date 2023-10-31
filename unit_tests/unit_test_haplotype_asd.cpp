//
// File: unit_test_haplotype_asd.cpp
// Date: 30 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit test the ASD calculation on haplotypes.
//

// Defines a genotype and ASD calculations.
#include "../src/Genotype.hpp"

// Used for testing operations.
//  This program has no output.
#include <cassert>

#include <iostream>
#include <bitset>
using namespace std;

int main()  {

    // Create two genotypes.
    //  This was already tested in unit_test_vcf_parser.cpp.
    Genotype a = GENOTYPE(1, 0);
    Genotype b = GENOTYPE(1, 2);
    Genotype g = GENOTYPE(1, 1);

    // Test GET_ALLELES.
    assert(GET_ALLELES(a) == 0x00000003);
    assert(GET_ALLELES(b) == 0x00000006);

    // Test ORIENTATION.
    assert(ORIENTATION(a) == 1);
    assert(ORIENTATION(b) == 0);

    // Test LSA.
    assert(LSA(a) == 0x00000001);
    assert(LSA(b) == 0x00000002);

    // Test LEFT_ALLELE.
    assert(LEFT_ALLELE(a) == 2);
    assert(LEFT_ALLELE(b) == 2);
    assert(LEFT_ALLELE(g) == 2);

    // Test RIGHT_ALLELE.
    assert(RIGHT_ALLELE(a) == 1);
    assert(RIGHT_ALLELE(b) == 4);
    assert(RIGHT_ALLELE(g) == 2);

    // Test UNORIENTED_ASD.
    assert(UNORIENTED_ASD(a, b) == 1);

    // Test ORIENTED_ASD.
    assert(ORIENTED_ASD(a, b) == 1);
    Genotype c = GENOTYPE(1, 3);
    Genotype f = GENOTYPE(1, 2);
    Genotype e = GENOTYPE(2, 2);
    assert(ORIENTED_ASD(b, c) == 1);
    assert(ORIENTED_ASD(b, f) == 0);
    assert(ORIENTED_ASD(a, e) == 2);

    // Create out samples' haplotypes.
    Genotype sample1[4];
    sample1[0] = GENOTYPE(1, 2);
    sample1[1] = GENOTYPE(1, 1);
    sample1[2] = GENOTYPE(2, 1);
    sample1[3] = GENOTYPE(2, 2);

    Genotype sample2[4];
    sample2[0] = GENOTYPE(1, 2);
    sample2[1] = GENOTYPE(1, 2);
    sample2[2] = GENOTYPE(2, 2);
    sample2[3] = GENOTYPE(2, 2);

    // Distance between samples starts off as 0.
    int d = 0;

    // Read in the first locus.
    assert((d = HAPLOTYPE_ASD(d, ORIENTED_ASD(sample1[0], sample2[0]), sample1[0], sample2[1])) == 0);
    
    // Read in the second locus.
    assert((d = HAPLOTYPE_ASD(d, ORIENTED_ASD(sample1[1], sample2[1]), sample1[1], sample2[1])) == 1);

    // Read in the third locus.
    assert((d = HAPLOTYPE_ASD(d, ORIENTED_ASD(sample1[2], sample2[2]), sample1[2], sample2[2])) == 2);

    // Read in the fourth locus.
    assert((d = HAPLOTYPE_ASD(d, ORIENTED_ASD(sample1[3], sample2[3]), sample1[3], sample2[3])) == 2);

}