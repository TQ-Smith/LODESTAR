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

int HAPLOTYPE_ASD(int asd, int left_a, int right_a, int left_b, int right_b) {
    if (asd == 2) {
        return 2;
    } else if (!(left_a ^ right_a ^ left_b ^ right_b)) {
        return 0;
    } else if (left_a != left_b && left_a != right_b && right_a != left_b && right_a != right_b ) {
        return 2;
    } else {
        return 1;
    }
}

int main()  {

    // Create two genotypes.
    //  This was already tested in unit_test_vcf_parser.cpp.
    Genotype a = GENOTYPE(0, 0);
    Genotype b = GENOTYPE(1, 1);
    Genotype c = GENOTYPE(0, 1);

    // Test genotype creation.
    assert(a == 0x80000001);
    assert(b == 0x80000002);
    assert(c == 0x00000003);

    // Test get alleles.
    assert(GET_ALLELES(a) == 0x00000001);
    assert(GET_ALLELES(b) == 0x00000002);
    assert(GET_ALLELES(c) == 0x00000003);

    // Test orientation.
    assert(ORIENTATION(a) == 1);
    assert(ORIENTATION(b) == 1);
    assert(ORIENTATION(c) == 0);

    // Test LSA.
    assert(LSA(a) == 0x00000001);
    assert(LSA(b) == 0x00000002);
    assert(LSA(c) == 0x00000001);

    // Test left allele.
    assert(LEFT_ALLELE(a) == 0x00000001);
    assert(LEFT_ALLELE(b) == 0x00000002);
    assert(LEFT_ALLELE(c) == 0x00000001);

    // Test right allele.
    assert(RIGHT_ALLELE(a) == 0x00000001);
    assert(RIGHT_ALLELE(b) == 0x00000002);
    assert(RIGHT_ALLELE(c) == 0x00000002);

    // Test unoriented ASD.
    assert(UNORIENTED_ASD(a, b) == 2);
    assert(UNORIENTED_ASD(a, c) == 1);
    assert(UNORIENTED_ASD(b, c) == 1);

    // Test oriented ASD.
    assert(ORIENTED_ASD(a, b) == 2);
    assert(ORIENTED_ASD(a, c) == 1);
    assert(ORIENTED_ASD(b, c) == 1);

    // Form part of the haplotype.
    Genotype a2 = GENOTYPE(1, 0);
    Genotype b2 = GENOTYPE(2, 1);
    Genotype c2 = GENOTYPE(1, 1);

    // Test left haplotype.
    assert(LEFT_HAPLOTYPE(a, a2) == 0x80000003);
    assert(LEFT_HAPLOTYPE(b, b2) == 0x80000006);
    assert(LEFT_HAPLOTYPE(c, c2) == 0x80000003);

    // Test right haplotype.
    assert(RIGHT_HAPLOTYPE(a, a2) == 0x80000001);
    assert(RIGHT_HAPLOTYPE(b, b2) == 0x80000002);
    assert(RIGHT_HAPLOTYPE(c, c2) == 0x80000002);

    // Test our haplotype asd calculation.
    int asd_ab = 0, asd_bc = 0, asd_ac = 0;
    assert(HAPLOTYPE_ASD(asd_ab, LEFT_HAPLOTYPE(a, a2), RIGHT_HAPLOTYPE(a, a2), LEFT_HAPLOTYPE(b, b2), RIGHT_HAPLOTYPE(b, b2)) == 2);
    assert(HAPLOTYPE_ASD(asd_bc, LEFT_HAPLOTYPE(b, b2), RIGHT_HAPLOTYPE(b, b2), LEFT_HAPLOTYPE(c, c2), RIGHT_HAPLOTYPE(c, c2)) == 1);
    assert(HAPLOTYPE_ASD(asd_ac, LEFT_HAPLOTYPE(a, a2), RIGHT_HAPLOTYPE(a, a2), LEFT_HAPLOTYPE(c, c2), RIGHT_HAPLOTYPE(c, c2)) == 1);

}