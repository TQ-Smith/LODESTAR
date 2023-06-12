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

    lodestar_pipeline("../50_of_CEU_chr2.vcf", "snp", 200, 200, 2);

}
