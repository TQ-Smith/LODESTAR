//
// File: VCFParser.cpp
// Date: 5 September 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines operations to parser a VCF file.
//

// Import the header.
#include "VCFParser.hpp"

// Used for cerr.
#include <iostream>

// Used for stringstream.
#include <sstream>

// A macro to hide the implementation of
//  creating a genotype.
#define GENOTPYE(a, b) (0x00000000 | 1 << a | 1 << b)

VCFParser::VCFParser(string file_name) {

    // Try to open file.
    // If file does not exist, print error and return.
    in_file = new ifstream;
    in_file -> open(file_name);
    if (!in_file -> is_open()) {
        cerr << "File " << file_name << " does not exist!" << endl;
        return;
    }

    // String buffer to hold line of the file.
    string buffer;

    // String to hold a sample's name.
    string sample;

    // Waste all of the meta-data at the beginning of the VCF file.
    //  When the header line is encountered, stop.
    while (buffer.substr(0, 2) != "#C") {
        getline(*in_file, buffer);
    }

    // Now, we iterate through the header line.
    //  Once the first 9 fields are encountered,
    //  we begin parsing the sample names.
    istringstream header(buffer);
    for (int i = 0; i < 9; i++) {
        getline(header, sample, '\t');
    }
    while (!header.eof()) {
        getline(header, sample, '\t');
        sample_names.push_back(sample);
    }

    // Finally, we set the number of samples.
    num_samples = sample_names.size();

}

bool VCFParser::getNextLocus(string* chromosome, int* position, bool* isMonomorphic, bool* isComplete, Genotype* genotypes) {

    // Set flags to true by default.
    *isMonomorphic = true;
    *isComplete = true;

    // A string buffer to read in a line.
    string buffer;

    // A string to hold the genotype of a sample.
    string sample;

    getline(*in_file, buffer);

    // If end of file reached, return false.
    if (in_file -> eof()) {
        return false;
    }

    // We iterate through the buffer and keep track of the current field.
    //  This is faster than creating a stringstream for each line.
    
    // The previous index in the buffer.
    int prev = -1;
    // The number of tabs encountered.
    int count = 0;

    // Used for genotype parsing.
    int delim;

    // Genotype to hold the previous sample's genotype.
    //  Used to keep track of monomorphic genotypes.
    Genotype previous_genotype;

    for (int i = 0; i <= (int) buffer.length(); i++) {
        if (i == (int) buffer.length() || buffer.at(i) == '\t') {

            // The first field is the chromosome.
            if (count == 0) {
                *chromosome = buffer.substr(prev + 1, i - prev - 1);
            }

            // The second field is the position.
            if (count == 1) {
                *position = stoi(buffer.substr(prev + 1, i - prev - 1));
            } 

            // The 9th field and greater hold the genotypes.
            if (count > 8) {

                // Get the genotype.
                //  NOTE: We are assuming the genotypes are first in the entry.
                sample = buffer.substr(prev + 1, i - prev);

                // If we encounter an incomplete genotype, the entry is incomplete,
                //  and we set the flag and exit the parsing. Includes case of just
                //  one allele.
                if (sample.length() == 2 || sample.at(0) == '.' || sample.at(2) == '.') {
                    *isComplete = false;
                    // Not EOF.
                    return true;
                }

                // Now, we parse the genotype.
                //  This should be a function but it will be called a lot
                //  and copy the string each time will be slow.

                // Get only the genotype information, which is assumed to be first.
                sample = sample.substr(0, sample.find(':'));

                // Find the allele seperator.
                if ((delim = sample.find('/')) == -1) {
                    delim = sample.find('|');
                }

                // Create and store our genotype.
                genotypes[count - 9] = GENOTPYE(stoi(sample.substr(0, delim)), stoi(sample.substr(delim + 1, sample.length() - delim)));

                // Test if monomorphic and keep track of previous genotype.
                if (count > 9 && genotypes[count - 9] != previous_genotype) {
                    *isMonomorphic = false;
                }
                previous_genotype = genotypes[count - 9];

            }

            prev = i;
            count++;
        }
    }

    // Not EOF.
    return true;

}

int VCFParser::getNumberOfSamples() {
    // Return the number of names.
    return sample_names.size();
}

string* VCFParser::getSampleNames() {
    // Create array to hold names.
    string* names = new string[num_samples];
    // Iterate through list and put each name in array.
    int i = 0;
    for (list<string>::iterator it = sample_names.begin(); it != sample_names.end(); it++) {
        names[i] = *it;
        i++;
    }
    return names;
}

VCFParser::~VCFParser() {
    // Close of in_file.
    in_file -> close();
    delete in_file;
}