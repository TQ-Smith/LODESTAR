
#include "RegionFilter.h"

#include <stdlib.h>

void insert_region(RegionFilter* filter, char* chrom, Region* region) {
    khiter_t k;
    int ret;
    k = kh_get(region, filter -> regions, chrom);
    if (k == kh_end(filter -> regions)) {
        k = kh_put(region, filter -> regions, chrom, &ret);
        kh_value(filter -> regions, k) = region;
    } else {
        
    }
}

bool parse_region(RegionFilter* filter, kstring_t* region, int startIndex, int endIndex) {
    char* leftEndpoint = NULL;
    char* rightEndpoint = NULL;
    int colonIndex = 0, hyphenIndex = 0, startLocus, endLocus;
    for (int i = startIndex; i <= endIndex; i++) {
        switch (ks_str(region)[i]) {
            case ':': colonIndex = i; break;
            case '-': hyphenIndex = i; break;
            default: break;
        }
    }
    startLocus = strtol(ks_str(region) + colonIndex + 1, &leftEndpoint, 10);
    endLocus = strtol(ks_str(region) + hyphenIndex + 1, &rightEndpoint, 10);

    Region* newRegion = (Region*) calloc(1, sizeof(Region));

    if (colonIndex == 0 && hyphenIndex == 0) {
        newRegion -> startLocus = 0;
        newRegion -> endLocus = ~0x0;
    } else if (hyphenIndex == 0 && leftEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = startLocus;
    } else if (hyphenIndex == endIndex && leftEndpoint == ks_str(region) + hyphenIndex) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = ~0x0;
    } else if (colonIndex + 1 == hyphenIndex && rightEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = 0;
        newRegion -> endLocus = -endLocus;
    } else if (leftEndpoint == ks_str(region) + hyphenIndex && rightEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = endLocus;
    } else {
        free(newRegion);
        return false;
    }

    char chrom[startIndex - endIndex + 1];
    for (int i = startIndex; i != colonIndex && i != endIndex + 1; i++)
        chrom[i - startIndex] = ks_str(region)[i];
    
    insert_region(filter, chrom, newRegion);

    return true;
}

bool parse_regions(RegionFilter* filter, kstring_t* inputRegions) {
    int prevIndex = -1;
    for (int i = 0; i <= ks_len(inputRegions); i++) {
        if (i == ks_len(inputRegions) || ks_str(inputRegions)[i] == ',') {
            if (!parse_region(filter, inputRegions, prevIndex + 1, i - 1))
                return false;
            prevIndex = i;
        }
    }
    return true;
}

RegionFilter* create_region_filter(kstring_t* inputRegions, bool takeComplement) {
    RegionFilter* filter = (RegionFilter*) calloc(1, sizeof(RegionFilter));
    filter -> takeComplement = takeComplement;
    filter -> regions = kh_init(region);
    if (parse_regions(filter, inputRegions)) {
        return filter;
    } else {
        destroy_region_filter(filter);
        return NULL;
    }
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
    kstring_t s;
    kputs("chr1:100-200,chr2,chr3:-200,chr4:300-,chr5:400,chr5-500", &s);
    create_region_filter(&s, false);
}