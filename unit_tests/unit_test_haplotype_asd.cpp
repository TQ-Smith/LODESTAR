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

/* Used to test the macro before I wrote the bit-op form.
int h_asd(int asd, int a1, int b1, int a2, int b2) {

    if (asd == 2) {
        return 2;
    } else if ((UNORIENTED_ASD(a1, b1) == 0) && (ORIENTATION(a1) != ORIENTATION(b1))) {
        return UNORIENTED_ASD(a2, b2);
    } else if (asd == 1) {
        return 2;
    } else {
        return ORIENTED_ASD(a2, b2);
    }

}
*/

int main()  {

    // Create two genotypes.
    //  This was already tested in unit_test_vcf_parser.cpp.
    Genotype a = GENOTYPE(0, 1);
    Genotype b = GENOTYPE(1, 0);
    Genotype c = GENOTYPE(1, 1);

    // Test get alleles.
    assert(GET_ALLELES(a) == 0x00000003);
    assert(GET_ALLELES(b) == 0x00000003);
    assert(GET_ALLELES(c) == 0x00000002);

    // Test orientation.
    assert(ORIENTATION(a) == 0);
    assert(ORIENTATION(b) == 1);
    assert(ORIENTATION(c) == 1);

    // Test LSA.
    assert(LSA(a) == 0x00000001);
    assert(LSA(b) == 0x00000001);
    assert(LSA(c) == 0x00000002);

    // Test left allele.
    assert(LEFT_ALLELE(a) == 0x00000001);
    assert(LEFT_ALLELE(b) == 0x00000002);
    assert(LEFT_ALLELE(c) == 0x00000002);

    // Test right allele.
    assert(RIGHT_ALLELE(a) == 0x00000002);
    assert(RIGHT_ALLELE(b) == 0x00000001);
    assert(RIGHT_ALLELE(c) == 0x00000002);

    // Test unoriented ASD.
    assert(UNORIENTED_ASD(a, b) == 0);
    assert(UNORIENTED_ASD(a, c) == 1);
    assert(UNORIENTED_ASD(b, c) == 1);

    // Test oriented ASD.
    assert(ORIENTED_ASD(a, b) == 2);
    assert(ORIENTED_ASD(a, c) == 1);
    assert(ORIENTED_ASD(b, c) == 1);

    // Create next loci in the haplotype.
    Genotype a2 = GENOTYPE(3, 3);
    Genotype b2 = GENOTYPE(3, 2);
    Genotype c2 = GENOTYPE(2, 2);

    // The asd values
    int ab_asd = 0, bc_asd = 0, ac_asd = 0;

    // Test haplotype asd. I *think* I got all cases.
    assert((ab_asd = HAPLOTYPE_ASD(ab_asd, a, b, a, b)) == 0);
    assert((bc_asd = HAPLOTYPE_ASD(bc_asd, b, c, b, c)) == 1);
    assert((ac_asd = HAPLOTYPE_ASD(ac_asd, a, c, a, c)) == 1);

    assert((ab_asd = HAPLOTYPE_ASD(ab_asd, a, b, a2, b2)) == 1);
    assert((bc_asd = HAPLOTYPE_ASD(bc_asd, b, c, b2, c2)) == 2);
    assert((ac_asd = HAPLOTYPE_ASD(ac_asd, a, c, a2, c2)) == 2);

    Genotype e = GENOTYPE(2, 3);
    Genotype g = GENOTYPE(2, 3);
    Genotype e2 = GENOTYPE(4, 5);
    Genotype g2 = GENOTYPE(5, 4);
    int eg_asd = 0;

    assert((eg_asd = HAPLOTYPE_ASD(eg_asd, e, g, e, g)) == 0);
    assert((eg_asd = HAPLOTYPE_ASD(eg_asd, e, g, e2, g2)) == 2);

    Genotype y = GENOTYPE(0, 1);
    Genotype z = GENOTYPE(1, 0);
    Genotype y2 = GENOTYPE(1, 0);
    Genotype z2 = GENOTYPE(0, 1);
    int yz_asd = 0;

    assert((yz_asd = HAPLOTYPE_ASD(yz_asd, y, z, y, z)) == 0);
    assert((yz_asd = HAPLOTYPE_ASD(yz_asd, y, z, y2, z2)) == 0);

}