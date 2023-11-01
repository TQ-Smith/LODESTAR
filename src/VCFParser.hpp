//
// File: VCFParser.hpp
// Date: 5 September 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines class to parser a VCF file.
//

#ifndef _VCF_PARSER_HPP_
#define _VCF_PARSER_HPP_

// Used to create genotypes.
#include "Genotype.hpp"

// Used to read in a compressed/uncompressed file.
#include "../lib/gzstream.h"

// We treat the lines of a file as a string.
#include <string>

// Our string operations lie in the std namespace.
using namespace std;

// Used to hold list of sample names.
#include <list>

class VCFParser {

    public:

        // Creates the parser object.
        //  Throws out all meta-information in VCF.
        //  Parses sample names and saves in list.
        //  Reader starts after header line.
        // Accepts:
        //  string file_name -> The name of the file to open.
        // Returns: void.
        VCFParser(string file_name);

        // Reads in next locus of VCF file.
        // Accepts:
        //  string* chromosome -> Sets the chromosome of the locus.
        //  int* position -> Sets the position of the locus.
        //  bool* isMonomorphic -> Sets boolean if locus was monomorphic across samples.
        //  bool* isComplete -> Sets if all samples had complete genotypes. Note, genotypes
        //                          are only parsed if all genotypes are complete.
        //  Genotype* genotypes -> Holds read in genotypes of samples at locus.
        // Returns: bool, True if not EOF, false otherwise.
        bool getNextLocus(string* chromosome, int* position, bool* isMonomorphic, bool* isComplete, Genotype* genotypes);

        // Get the number of samples.
        // Accepts: void.
        // Returns: int, The number of samples.
        int getNumberOfSamples();

        // Returns an array of the sample names.
        // Accepts: void.
        // Returns: string*, The array of the sample names.
        string* getSampleNames();

        // Tests if file is open.
        // Accepts: void.
        // Returns: bool, If the file is open.
        bool isOpen();

        // Destroy the parser object.
        // Accepts: void.
        // Returns: void.
        ~VCFParser();

    private:

        // The name of the VCF file.
        string file_name;

        // The number of samples in the file.
        int num_samples;

        // List of the sample names.
        list<string> sample_names;

        // Our in file stream object to read in the contents of the file.
        igzstream* in_file;

};

#endif