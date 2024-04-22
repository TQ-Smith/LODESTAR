
#ifndef _REGION_FILTER_H_
#define _REGION_FILTER_H_

#include "../lib/kstring.h"

#include "../lib/khash.h"

#include "stdbool.h"

typedef struct region {
    unsigned int startLocus;
    unsigned int endLocus;
    struct region* next;
} Region;

KHASH_MAP_INIT_STR(region, Region*)

typedef struct {
    bool takeComplement;
    khash_t(region)* regions;
} RegionFilter;

RegionFilter* create_region_filter(kstring_t* inputRegions, bool takeComplement);

bool query_locus(RegionFilter* filter, kstring_t* chrom, unsigned int locus);

bool query_region(RegionFilter* filter, kstring_t* chrom, unsigned int startLocus, unsigned int endLocus);

void destroy_region_filter(RegionFilter* filter);

#endif