//
// File: Genotype.hpp
// Date: 30 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Define implementation of representing genotypes.
//

#ifndef _GENOTYPE_HPP_
#define _GENOTYPE_HPP_

// We typedef int to represent a genotype at a locus.
//  This slight distinction should help readability.
typedef int Genotype;

// A macro to hide the implementation of creating a genotype. 
//  Sets bits corresponding to the allele number. Max of 31 alleles = sizeof(int) - 1.
//  Sets MSB if the greater allele is the left allele. Helps maintain
//  order when creating haplotypes. So, 1/0 would have orientation set to 1.
#define GENOTYPE(a, b) (0x00000000 | 1 << a | 1 << b | ((a >= b) << 31))

// A macro to get just the alleles without orientation.
#define GET_ALLELES(a) (a & 0x7FFFFFFF)

// A macro to get the orientation of a genotype.
#define ORIENTATION(a) ((a & 0x80000000) >> 31)

// A macro to get an integer with just the right-most set bit.
//  Instead of least-significant-bit, we call it the least-significant-allele.
#define LSA(a) (a - (a & a - 1))

// A macro to get the oriented left allele. For example, if the genotype is 2/0 
//  with orientation, the macro returns the number 0x00000004.
#define LEFT_ALLELE(a) (ORIENTATION(a) * (!(GET_ALLELES(a) ^ LSA(a)) * GET_ALLELES(a) + !!(GET_ALLELES(a) ^ LSA(a)) * (GET_ALLELES(a) - LSA(a))) + !ORIENTATION(a) * LSA(a))

// A macro to get the oriented right allele. For example, if the genotype is 2/0 
//  with orientation, the macro returns the number 0x00000000.
#define RIGHT_ALLELE(a) (!ORIENTATION(a) * (GET_ALLELES(a) - LSA(a)) + ORIENTATION(a) * LSA(a))

// A macro to calculate the ASD between two genotypes IGNORING orientation.
//  ASD in the sense of the number of different alleles between the two genotypes.
#define UNORIENTED_ASD(a, b) (!!(GET_ALLELES(a) ^ GET_ALLELES(b)) + !(GET_ALLELES(a) & GET_ALLELES(b)))

// A macro to calculate the number of shared allele between two genotypes WITH orientation.
#define ORIENTED_ASD(a, b) (!!(LEFT_ALLELE(a) ^ LEFT_ALLELE(b)) + !!((RIGHT_ALLELE(a) ^ RIGHT_ALLELE(b))))

// A macro to form the left haplotype between two genotype.
//  If the greater of the two alleles comes last, set most significant bit.
#define LEFT_HAPLOTYPE(a, b) (LEFT_ALLELE(a) | LEFT_ALLELE(b) | ((LEFT_ALLELE(b) >= LEFT_ALLELE(a)) << 31))

// A macro to form the right haplotype between two genotype.
//  If the greater of the two alleles comes last, set most significant bit.
#define RIGHT_HAPLOTYPE(a, b) (RIGHT_ALLELE(a) | RIGHT_ALLELE(b) | ((RIGHT_ALLELE(b) >= RIGHT_ALLELE(a)) << 31))

// We define out haplotype asd function.
//  NOTE: This only works for two loci.
//        Counts the number of shared haplotypes.
//  Equivalent to the following if-statement:
//  
//  if (left_a == left_b && right_a == right_b)
//      return 0;
//  else if ((left_a != left_b) && (left_a != left_b) && (right_a != left_b) && (right_a != right_b))
//      return 2;
//  else
//      return 1;
//
// Accepts:
//  left_a -> The left haplotype of sample a.
//  right_a -> The right haplotype of sample a.
//  left_b -> The left haplotype of sample b.
//  right_b -> The right haplotype of sample b.
// Returns: double, The ASD between the two haplotypes.
#define HAPLOTYPE_ASD(left_a, right_a, left_b, right_b) (\
    (!(left_a ^ left_b) && !(right_a ^ right_b)) * 0.0 +\
    ((!!(left_a ^ left_b) || !!(right_a ^ right_b)) && (!!(left_a ^ left_b) && !!(left_a ^ left_b) && (!!right_a ^ left_b) && !!(right_a ^ right_b))) * 2.0 + \
    ((!!(left_a ^ left_b) || !!(right_a ^ right_b)) && (!(left_a ^ left_b) || !(left_a ^ left_b) || !(right_a ^ left_b) || !(right_a ^ right_b))) * 1.0 \
)

#endif