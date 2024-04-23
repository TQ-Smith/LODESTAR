
#include "RegionFilter.h"

void insert_region(RegionFilter* filter, const char* chrom, Region* region) {
    khint_t k;
    int absent;
    k = kh_get(region, filter -> regions, chrom); 
    if (k == kh_end(filter -> regions)) {
        // printf("Insert %s: %u to %u\n", chrom, region -> startLocus, region -> endLocus);
        k = kh_put(region, filter -> regions, chrom, &absent);
        if (absent) 
            kh_key(filter -> regions, k) = strdup(chrom);
        kh_value(filter -> regions, k) = region;
    } else {
        // printf("Append %s: %u to %u\n", chrom, region -> startLocus, region -> endLocus);
        Region* current = kh_value(filter -> regions, k);
        Region* previous = current;
        while (current != NULL && region -> startLocus > current -> startLocus) {
            previous = current;
            current = current -> next;
        }
        // Note this is imprefect, since when much larger intervals are added, the subset intervals are not purged.
        //  This should not matter though.
        if (current != NULL && region -> endLocus <= current -> endLocus) {
            free(region);
        } else {
            previous -> next = region;
            region -> next = current;
        }
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
    newRegion -> next = NULL;

    if (colonIndex == 0 && hyphenIndex == 0) {
        newRegion -> startLocus = 0;
        newRegion -> endLocus = (unsigned int) ~0x0;
    } else if (hyphenIndex == 0 && leftEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = startLocus;
    } else if (hyphenIndex == endIndex && leftEndpoint == ks_str(region) + hyphenIndex) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = (unsigned int) ~0x0;
    } else if (colonIndex + 1 == hyphenIndex && rightEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = 0;
        newRegion -> endLocus = endLocus;
    } else if (leftEndpoint == ks_str(region) + hyphenIndex && rightEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = endLocus;
    } else {
        free(newRegion);
        return false;
    }

    if (newRegion -> startLocus > newRegion -> endLocus) {
        free(newRegion);
        return false;
    }

    int chromEndIndex = colonIndex != 0 ? colonIndex - 1 : endIndex;
    char chrom[chromEndIndex - startIndex + 2];
    for (int i = 0; i < chromEndIndex - startIndex + 1; i++)
        chrom[i] = ks_str(region)[startIndex + i];
    chrom[chromEndIndex - startIndex + 1] = '\0';

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
    if (parse_regions(filter, inputRegions))
        return filter;
    destroy_region_filter(filter);
    return NULL;
}

bool query_locus(RegionFilter* filter, kstring_t* chrom, unsigned int locus) {
    khint_t k = kh_get(region, filter -> regions, ks_str(chrom));
    if (k != kh_end(filter -> regions)) {
        Region* current = kh_value(filter -> regions, k);
        while (current != NULL && locus >= current -> startLocus) {
            if (locus <= current -> endLocus)
                return filter -> takeComplement ^ true;
            current = current -> next;
        }
    }
    return filter -> takeComplement ^ false;
}

bool query_overlap(RegionFilter* filter, kstring_t* chrom, unsigned int startLocus, unsigned int endLocus) {
    khint_t k = kh_get(region, filter -> regions, ks_str(chrom));
    if (k != kh_end(filter -> regions)) {
        Region* leftRegion = kh_value(filter -> regions, k);
        while (leftRegion != NULL && startLocus >= leftRegion -> startLocus)
            leftRegion = leftRegion -> next;

        Region* rightRegion = leftRegion;
        while (rightRegion != NULL && endLocus >= rightRegion -> endLocus)
            rightRegion = rightRegion -> next;

        if (leftRegion == NULL || (startLocus > leftRegion -> endLocus && endLocus < rightRegion -> startLocus))
            return filter -> takeComplement ^ false;
        else 
            return filter -> takeComplement ^ true;
    }
    return filter -> takeComplement ^ false;
}

void destroy_region_filter(RegionFilter* filter) {
    if (filter == NULL)
        return;
    Region* current = NULL;
    Region* temp = NULL;
    khint_t j;
    for (j = kh_begin(filter -> regions); j != kh_end(filter -> regions); j++) {
		if (!kh_exist(filter -> regions, j))
            continue;
		current = kh_val(filter -> regions, j);
        while (current != NULL) {
            // printf("Free %s: %u to %u\n", kh_key(filter -> regions, j), current -> startLocus, current -> endLocus);
            temp = current;
            current = current -> next;
            free(temp);
        }
        free((char*) kh_key(filter -> regions, j));
	}
    kh_destroy(region, filter -> regions);
    free(filter);
}

/*
int main() {
    kstring_t* s = (kstring_t*) calloc(1, sizeof(kstring_t));
    kstring_t* s1 = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs("chr3", s1);
    kstring_t* s2 = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs("chr4", s2);
    kstring_t* s3 = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs("chr5", s3);
    kstring_t* s4 = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs("chr1", s4);
    printf("Create Regions:\n");
    kputs("chr1:100-200,chr2,chr3:-200,chr4:300-,chr5:400,chr5:500-600,chr1:1000-2000,chr1:500-600,chr2:100-200", s);
    RegionFilter* filter = create_region_filter(s, false);
    printf("\n");
    printf("Query:\n");
    printf("Query locus chr3:100: %d\n", query_locus(filter, s1, 100));
    printf("Query locus chr4:100: %d\n", query_locus(filter, s2, 100));
    printf("Query locus chr5:550: %d\n", query_locus(filter, s3, 550));
    printf("Query overlap chr1:900-2200: %d\n", query_overlap(filter, s4, 900, 2200));
    printf("Query overlap chr1:150-250: %d\n", query_overlap(filter, s4, 150, 250));
    printf("Query overlap chr1:50-150: %d\n", query_overlap(filter, s4, 50, 150));
    printf("Query overlap chr1:3000-4000: %d\n", query_overlap(filter, s4, 3000, 4000));
    printf("\n");
    printf("Free Regions\n");
    free((char*) ks_release(s)); free(s);
    free((char*) ks_release(s1)); free(s1);
    free((char*) ks_release(s2)); free(s2);
    free((char*) ks_release(s3)); free(s3);
    free((char*) ks_release(s4)); free(s4);
    destroy_region_filter(filter);
}
*/