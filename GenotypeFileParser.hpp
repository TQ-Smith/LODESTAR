//
// File: GenotypeFileParser.hpp
// Started: 22 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Header file for GenotypeFileParser.cpp
//

#ifndef GENOTYPE_FILE_PARSER_HPP_
#define GENOTYPE_FILE_PARSER_HPP_

// I try to keep other includes out of the header, but here, it is necessary.

// Used to read in a file.
#include <fstream>

// We treat the lines of a file as a string.
#include <string>

// Our string operations lie in the std namespace.
using namespace std;

struct Genotype {
    int allele;
};

// Opens a file.
// Accepts:
//  ifstream& in_file -> A preinitialized in-file stream.
//  string file_name -> The name of the file to open.
// Returns: void.
void open_file(ifstream& in_file, string file_name);

// Closes a file.
// Accepts:
//  ifstream& in_file -> A preinitalized in-file stream.
// Returns: void.
void close_file(ifstream& in_file);

// From an opened VCF file, get sample names and the number of samples.
// Accepts:
//  ifstream& in_file -> An in-file stream to a VCF file.
//  string*& sample_names -> A NULL string pointer. Sets pointer to a new string,
//                              where each element is the name of a sample.
//  int& n -> Sets integer to the number of sample.
// Returns: void.
void get_sample_names(ifstream& in_file, string*& sample_names, int& n);

// Get the next locus from the supplied VCF file.
//  NOTE: Assumes get_sample_names have already been run.
//          That is, the stream is pointing to the first non-meta line of the VCF file.
// Accepts:
//  ifstream& in_file -> An in-file stream to a VCF file.
//  string& chrom -> Sets the name of the locus's chromosome.
//  int& position -> Sets the locus's position.
//  Genotype* genotypes -> Sets the genotype for easch sample.
//  bool& isComplete -> Sets if every sample's genotype was complete.
//  bool& isEOF -> Sets if the next line in the stream was end of the file.
//  int n -> The number of sample.
void get_next_loci(ifstream& in_file, string& chrom, int& position, Genotype* genotypes, bool& isComplete, bool& isMonomorphic, bool& isEOF, int n);

#endif