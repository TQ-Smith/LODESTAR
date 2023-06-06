//
// File: lodestar.cpp
// Started: 17 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Parses command line arguments and acts as the main program.
//

#include "CommandLineArgumentParser.hpp"

#include "Engine.hpp"

int main() {

    lodestar_pipeline("/home/tqs5778/Documents/Data/50_of_CEU_CHB_YRI_chr2.vcf", "snp", 50, 50, 2);

}