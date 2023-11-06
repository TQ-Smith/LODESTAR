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

    // Test unoriented DISSIMILARITY.
    assert(UNORIENTED_DISSIMILARITY(a, b) == 2);
    assert(UNORIENTED_DISSIMILARITY(a, c) == 1);
    assert(UNORIENTED_DISSIMILARITY(b, c) == 1);

    // Test oriented DISSIMILARITY.
    assert(ORIENTED_DISSIMILARITY(a, b) == 2);
    assert(ORIENTED_DISSIMILARITY(a, c) == 1);
    assert(ORIENTED_DISSIMILARITY(b, c) == 1);

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
    //  NOTE: HAPLOTYPE_DISSIMILARITY only considers haplotypes of two loci.
    // Two different haplotypes.
    assert(HAPLOTYPE_DISSIMILARITY(LEFT_HAPLOTYPE(a, a2), RIGHT_HAPLOTYPE(a, a2), LEFT_HAPLOTYPE(b, b2), RIGHT_HAPLOTYPE(b, b2)) == 2);
    // Test only one in common.
    assert(HAPLOTYPE_DISSIMILARITY(LEFT_HAPLOTYPE(b, b2), RIGHT_HAPLOTYPE(b, b2), LEFT_HAPLOTYPE(c, c2), RIGHT_HAPLOTYPE(c, c2)) == 1);
    assert(HAPLOTYPE_DISSIMILARITY(LEFT_HAPLOTYPE(a, a2), RIGHT_HAPLOTYPE(a, a2), LEFT_HAPLOTYPE(c, c2), RIGHT_HAPLOTYPE(c, c2)) == 1);
    // Test two identical haplotypes.
    assert(HAPLOTYPE_DISSIMILARITY(LEFT_HAPLOTYPE(a, a2), RIGHT_HAPLOTYPE(a, a2), LEFT_HAPLOTYPE(a, a2), RIGHT_HAPLOTYPE(a, a2)) == 0);

}