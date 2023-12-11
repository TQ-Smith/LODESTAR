
#include "VCFGenotypeParser.h"

VCFGenotypeParser* init_vcf_genotype_parser(char* file_name) {

    BGZF* file = bgzf_open(file_name, 'r');

    if (bgzf_check_EOF(file)) {
        return NULL;
    }

}

bool getNextLocus(VCFGenotypeParser* parser, char* chromosome, int* position, GENOTYPE* genotypes) {
    
}

char** getSampleNames(VCFGenotypeParser* parser) {

}

destroy_vcf_genotype_parser(VCFGenotypeParser* parser) {

}

int main() {

}