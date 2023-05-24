//
// File: GenotypeFileParser.cpp
// Started: 22 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains methods to parse various file types.
//

// #include "GenotypeFileParser.hpp"
//

#include <iostream>
#include <fstream>
#include <string>
#include <list>
using namespace std;

void open_file(ifstream& in_file, string file_name) {
    in_file.open(file_name);
}

void get_sample_names(ifstream& in_file, string*& sample_names, int& n) {
    
    string line;

    list<string> samples;
    int count = 0;

    while (line.substr(0, 2) != "#C") {
        getline(in_file, line);
    }

    int prev = -1;
    for(int i = 0; i < line.size(); i++) {
        if (line.at(i) == '\t') {
            samples.push_back(line.substr(prev + 1, i - prev));
            prev = i;
            count++;
        }
    }
    
    
    n = count - 9;

    sample_names = new string[n];

    for (int i = 0; i < 9; i++) {
        samples.pop_front();
    }

    int k = 0;
    for (string i: samples) {
        sample_names[k] = i;
        k++;
    }

}

void close_file(ifstream& in_file) {
    in_file.close();
}

struct Genotype {
    int chr1;
    int chr2;
};

void get_next_loci(ifstream& in_file, string& chrom, int& position, bool& success, Genotype* genotypes, int n) {

    // Change this!!!!

    string line, sample;

    int colon;

    getline(in_file, line);

    cout << line << endl << endl;

    int length = line.length();

    int tab = line.find('\t');

    chrom = line.substr(0, tab);

    line = line.substr(tab + 1, length);

    tab = line.find('\t');

    position = stoi(line.substr(0, tab));

    for (int i = 0; i < 8; i++) {
        tab = line.find('\t');
        line = line.substr(tab + 1, length);
    }

    for (int i = 0; i < n; i++) {
        tab = line.find('\t');
        sample = line.substr(0, tab);

        if ((colon = sample.find(':')) != -1) {
            sample = sample.substr(0, colon);
        }

        if ((colon = sample.find('/')) == -1) {
            colon = sample.find('|');
        }

        
        string chr1 = sample.substr(0, colon);
        if (chr1 == ".") {
            success = false;
            break;
        }

        genotypes[i].chr1 = stoi(chr1);
        genotypes[i].chr2 = stoi(sample.substr(colon + 1, sample.length()));
        
        line = line.substr(tab + 1, length);
    }

}


int main() {
    ifstream in_file;
    open_file(in_file, "dingo.vcf");

    string* samples = NULL;
    int n;

    get_sample_names(in_file, samples, n);

    Genotype* sample_genotypes = new Genotype[n];

    string chrom;
    int position;
    bool success;

    get_next_loci(in_file, chrom, position, success, sample_genotypes, n);

    for (int i = 0; i < n; i++) {
        cout << sample_genotypes[i].chr1 << " " << sample_genotypes[1].chr2 << endl;
    }

    delete [] sample_genotypes;

    close_file(in_file);
}
