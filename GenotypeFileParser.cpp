//
// File: GenotypeFileParser.cpp
// Started: 22 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains methods to parse various file types.
//

#include "GenotypeFileParser.hpp"


#include <iostream>
using namespace std;

// We use lists read in the sample names.
#include <list>

void open_file(ifstream& in_file, string file_name) {
    // Right now, just open the file.
    //  Later, add error handeling.
    in_file.open(file_name);
}

void get_sample_names(ifstream& in_file, string*& sample_names, int& n) {
    
    // Buffer to read in a line of a file.
    string line;

    // A list to hold the sample names.
    list<string> samples;

    // A counter for the number of fields in the VCF file.
    int count = 0;

    // Waste all of the meta-data at the beginning of the VCF file.
    //  When the header line is encountered, stop.
    while (line.substr(0, 2) != "#C") {
        getline(in_file, line);
    }

    // Get the length of the line.
    int length = line.length();

    // We split the header line on the tab character.
    //  Once the first 9 tabs were encountered, we have reached the
    //  sample names. Save each sample name in the list.
    int prev = -1;
    for(int i = 0; i < length; i++) {
        if (line.at(i) == '\t') {
            if (count > 8) {
                samples.push_back(line.substr(prev + 1, i - prev));
            }
            prev = i;
            count++;
        }
    }

    // Add the last sample.
    samples.push_back(line.substr(prev + 1, length - prev));
    
    // Subtract the 9 non-sample related fields.
    n = count - 8;

    // Create our array for sample names.
    sample_names = new string[n];

    // Store the list of sample names into the arrray.
    int k = 0;
    for (string i: samples) {
        sample_names[k] = i;
        k++;
    }

}

// A helper function to parse a genotype for a sample.
// Accepts:
//  string sample -> The genotype to parse.
//  Genotype& locus -> Set the genotype.
// Returns: void.
void parse_genotype(string sample, Genotype& locus) {
    // A string to hold one of the haploid values.
    string haploid;

    // Used to remove ':'s in the VCF genotype field and used
    //  to split genotypes into integers.
    int delim_index;

    // Just get the genotype in the entry if ':' is present.
    if ((delim_index = sample.find(':')) != -1) {
        sample = sample.substr(0, delim_index);
    }

    // Split genotype by '/' or '|'.
    if ((delim_index = sample.find('/')) == -1) {
        delim_index = sample.find('|');
    }

    // Get the first haploid.
    haploid = sample.substr(0, delim_index);

    // Store it in the sample's genotype.
    locus.chr1 = stoi(haploid);

    // Get the second haploid and store it in the sample's genotype.
    locus.chr2 = stoi(sample.substr(delim_index + 1, sample.length() - delim_index));
}

// I really don't love how I wrote this method.
//  Keep back in mind to fix later.
//  This also needs to be fast.
//  If we used a buffer instead of getline, this would be faster.
void get_next_loci(ifstream& in_file, string& chrom, int& position, Genotype* genotypes, bool& isComplete, bool& isEOF, int n) {

    // When we start, it isn't EOF;
    isEOF = false;
    
    // Buffer to read in a line of the file.
    string line;

    // Read in the chromosome name.
    //  If EOF, set flag and exit method.
    if (!getline(in_file, line)) {
        isEOF = true;
        return;
    }

    // Get the length of the rest of the line.
    int length = line.length();

    // A string to hold a sample's genotype.
    string sample;

    // We split the line on '\t', but throw away the fields that
    //  do not correspond to sample.
    // NOTE: I don't think I can get around the deep nesting.
    int prev = -1;
    int count = 0;
    for (int i = 0; i < length; i++) {
        if (line.at(i) == '\t') {
            // First field is the chromosome.
            if (count == 0) {
                chrom = line.substr(prev + 1, i - prev);
            // Second field is the position.
            } else if (count == 1) {
                position = stoi(line.substr(prev + 1, i - prev));
            // Eight fields and greater are the sample's genotypes.
            } else if (count > 8) {

                // Get the genotype.
                //  NOTE: We are assuming the genotypes are first in the entry.
                sample = line.substr(prev + 1, i - prev);

                // If we encounter an incomplete genotype, the entry is incomplete,
                //  and we set the flag and exit the parsing.
                if (sample.at(0) == '.') {
                    isComplete = false;
                    return;
                }

                parse_genotype(sample, genotypes[count - 9]);
                
            }
            prev = i;
            count++;
        }
    }

    // Parse the last genotype.
    sample = line.substr(prev + 1, length - prev);

    // If we encounter an incomplete genotype, the entry is incomplete,
    //  and we set the flag and exit the parsing.
    if (sample.at(0) == '.') {
        isComplete = false;
        return;
    }
    parse_genotype(sample, genotypes[count - 9]);

    // If we go through the whole file, then the record was complete.
    isComplete = true;

}

void close_file(ifstream& in_file) {
    // Right now, just open the file.
    //  Later, add error handeling.
    in_file.close();
}