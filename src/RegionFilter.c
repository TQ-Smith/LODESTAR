
#include "RegionFilter.h"

#include <stdlib.h>

RegionFilter* create_region_filter(bool takeComplement) {
    RegionFilter* filter = (RegionFilter*) calloc(1, sizeof(RegionFilter));
    filter -> takeComplement = takeComplement;
    filter -> regions = kh_init(region);
    return filter;
}

Region* parse_region(kstring_t* region, int startIndex, int endIndex) {

}

int parse_regions(RegionFilter* filter, kstring_t* inputRegions) {
    int prevIndex = 0;
    for (int i = 0; i <= ks_len(inputRegions); i++) {
        if (i == ks_len(inputRegions) || ks_str(inputRegions)[i] == ',') {
            


            prevIndex = i;
        }
    }
    return 0;
}

bool query_locus(RegionFilter* filter, kstring_t* chrom, unsigned int locus) {

}

bool query_region(RegionFilter* filter, kstring_t* chrom, unsigned int startLocus, unsigned int endLocus) {

}

void destroy_region_filter(RegionFilter* filter) {
    if (filter == NULL)
        return;
    kh_destroy(region, filter -> regions);
    free(filter);
}

int main() {
    
}