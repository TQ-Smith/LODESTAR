
// File: RegionSet.c
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Used to define a set of genomic coordinates.

#include "RegionSet.h"

// Insert and interval on a chromosome into the set.
// Accepts:
//  RegionSet* set -> The set to insert the region.
//  const char* chrom -> The chromosome the interval lies on.
//  Region* region -> The region that represents an interval on a chromosome.
// Returns: void.
void insert_region(RegionSet_t* set, const char* chrom, Region_t* region) {
    khint_t k;
    int absent;
    k = kh_get(Region, set -> regions, chrom); 
    // If the chromosome is not in the hash table...
    if (k == kh_end(set -> regions)) {
        // Insert the chromosome as a key.
        k = kh_put(Region, set -> regions, chrom, &absent);
        // If the bucket is empty, copy the chromosome key into the bucket.
        if (absent) 
            kh_key(set -> regions, k) = strdup(chrom);
        // Insert the region.
        kh_value(set -> regions, k) = region;
    } else {
        Region_t* current = kh_value(set -> regions, k);
        Region_t* previous = current;
        // Iterate through the sorted linked-list to insert the new interval.
        while (current != NULL && region -> startLocus > current -> startLocus) {
            previous = current;
            current = current -> next;
        }
        // If the interval is contained in an already existant larger interval, no need to add the region.
        //  Note, we can make this more robust by merging overlapping intervals into a larger interval.
        if (current != NULL && region -> endLocus <= current -> endLocus) {
            free(region);
        } else {
            // Link region into the list.
            previous -> next = region;
            region -> next = current;
        }
    }
}

// Parse an individual region.
//  RegionSet_t* set -> The set to search.
//  kstring_t* region -> Our string that defines a single region. Follows grammar defined in RegionSet.h
//  int startIndex -> The index of the first character defining the region in the region string.
//  int endIndex -> The index of the last character defining the region in the region string.
// Returns: bool, True if region was syntactically valid. False otherwise.
bool parse_region(RegionSet_t* set, kstring_t* region, int startIndex, int endIndex) {
    char* leftEndpoint = NULL;
    char* rightEndpoint = NULL;
    int colonIndex = 0, hyphenIndex = 0, startLocus, endLocus;
    // Find the position of the colon and hyphen in the region, if they exist.
    for (int i = startIndex; i <= endIndex; i++) {
        switch (ks_str(region)[i]) {
            case ':': colonIndex = i; break;
            case '-': hyphenIndex = i; break;
            default: break;
        }
    }
    // startLocus is the start coordinate of the inveral.
    startLocus = strtol(ks_str(region) + colonIndex + 1, &leftEndpoint, 10);
    // startLocus is the end coordinate of the inveral.
    endLocus = strtol(ks_str(region) + hyphenIndex + 1, &rightEndpoint, 10);

    // Note, leftEndpoint and rightEndpoint point to the starting character after parsing either integer.

    // Create our region.
    Region_t* newRegion = (Region_t*) calloc(1, sizeof(Region_t));
    newRegion -> next = NULL;

    // If neither a colon or hyphen were encountered, the region defines a whole chromosome.
    if (colonIndex == 0 && hyphenIndex == 0) {
        newRegion -> startLocus = 0;
        newRegion -> endLocus = (unsigned int) ~0x0;
    // If not hyphen was encountered and the startLocus was successfully parsed, we are dealing with the STR:INT case.
    //  Region is just one locus.
    } else if (hyphenIndex == 0 && leftEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = startLocus;
    // If the hyphen preceeds the start of the next region, we are dealing with the case STR:INT-
    } else if (hyphenIndex == endIndex && leftEndpoint == ks_str(region) + hyphenIndex) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = (unsigned int) ~0x0;
    // If the hyphen comes after the colon, we are dealing with the case STR:-INT
    } else if (colonIndex + 1 == hyphenIndex && rightEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = 0;
        newRegion -> endLocus = endLocus;
    // If both the colon and hyphen are presenent with both the startLocus and endLocus, we are dealing with
    //  the case STR:INT-INT.
    } else if (leftEndpoint == ks_str(region) + hyphenIndex && rightEndpoint == ks_str(region) + endIndex + 1) {
        newRegion -> startLocus = startLocus;
        newRegion -> endLocus = endLocus;
    // If none of the above conditions are met, we free our region and return false, indicating the region
    //  was not syntatically valid.
    } else {
        free(newRegion);
        return false;
    }

    // If the start coordinate comes after the end coordinate, the region is not valid. 
    //  This is semantically incorrect. Return error.
    if (colonIndex != 0 && hyphenIndex != 0 && newRegion -> startLocus > newRegion -> endLocus) {
        free(newRegion);
        return false;
    }

    // Copy chromosome to a new string to be used as a key in the hash table.
    int chromEndIndex = colonIndex != 0 ? colonIndex - 1 : endIndex;
    char chrom[chromEndIndex - startIndex + 2];
    for (int i = 0; i < chromEndIndex - startIndex + 1; i++)
        chrom[i] = ks_str(region)[startIndex + i];
    chrom[chromEndIndex - startIndex + 1] = '\0';

    // Insert the region into the set.
    insert_region(set, chrom, newRegion);

    // Successfully parsed region.
    return true;
}

// Parse a list of REGION.
// Accepts:
//  RegionSet_t* set -> The set to search.
//  kstring_t* inputRegions -> Our string that defines intervals. Follows grammar defined in RegionSet.h
// Returns: bool, True if all regions were syntactically valid. False otherwise.
bool parse_regions(RegionSet_t* set, kstring_t* inputRegions) {
    int prevIndex = -1;
    for (int i = 0; i <= ks_len(inputRegions); i++) {
        if (i == ks_len(inputRegions) || ks_str(inputRegions)[i] == ',') {
            if (!parse_region(set, inputRegions, prevIndex + 1, i - 1))
                return false;
            prevIndex = i;
        }
    }
    return true;
}

RegionSet_t* init_region_set(kstring_t* inputRegions, bool takeComplement) {
    // Allocate a set.
    RegionSet_t* set = (RegionSet_t*) calloc(1, sizeof(RegionSet_t));
    set -> takeComplement = takeComplement;
    set -> regions = kh_init(Region);
    if (parse_regions(set, inputRegions))
        return set;
    // If we were unable to parse our inputRegion string, destroy our created set and return NULL.
    destroy_region_set(set);
    return NULL;
}

bool query_locus(RegionSet_t* set, kstring_t* chrom, unsigned int locus) {
    printf("%d\n", locus);
    khint_t k = kh_get(Region, set -> regions, ks_str(chrom));
    // If the chromosome is in the hash set ...
    if (k != kh_end(set -> regions)) {
        // Iterate through intervals on the chromosome.
        Region_t* current = kh_value(set -> regions, k);
        while (current != NULL && locus >= current -> startLocus) {
            // Coordinate is in the interval 
            if (locus <= current -> endLocus)
                return set -> takeComplement ^ true;
            current = current -> next;
        }
    }
    // Chromosome is not in the hash table.
    return set -> takeComplement ^ false;
}

bool query_overlap(RegionSet_t* set, kstring_t* chrom, unsigned int startLocus, unsigned int endLocus) {
    khint_t k = kh_get(Region, set -> regions, ks_str(chrom));
    // If the chromosome is in the hash set ...
    if (k != kh_end(set -> regions)) {
        Region_t* leftRegion = NULL;
        Region_t* rightRegion = NULL;
        // Search if startLocus lands in an interval.
        leftRegion = kh_value(set -> regions, k);
        while (leftRegion != NULL && startLocus >= leftRegion -> startLocus) {
            if (startLocus <= leftRegion -> endLocus)
                return set -> takeComplement ^ true;
            leftRegion = leftRegion -> next;
        }
        // Search if endLocus lands in an interval OR a whole interval is contained in the window.
        rightRegion = leftRegion;
        while (rightRegion != NULL && endLocus >= rightRegion -> startLocus) {
            if (startLocus <= rightRegion -> endLocus || endLocus <= rightRegion -> endLocus)
                return set -> takeComplement ^ true;
            rightRegion = rightRegion -> next;
        }
        
    }
    // Chromosome is not in the hash table.
    return set -> takeComplement ^ false;
}

void destroy_region_set(RegionSet_t* set) {
    if (set == NULL)
        return;
    Region_t* current = NULL;
    Region_t* temp = NULL;
    khint_t j;
    // Iterate through each chromosome.
    for (j = kh_begin(set -> regions); j != kh_end(set -> regions); j++) {
		if (!kh_exist(set -> regions, j))
            continue;
        // If bucket is not empty, then destroy the sorted linked-list.
		current = kh_val(set -> regions, j);
        while (current != NULL) {
            temp = current;
            current = current -> next;
            free(temp);
        }
        // We must destroy the keys since they are strings.
        free((char*) kh_key(set -> regions, j));
	}
    // Destroy the hash table and the set itself.
    kh_destroy(Region, set -> regions);
    free(set);
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
    RegionSet_t* set = init_region_set(s, false);
    printf("\n");
    printf("Query:\n");
    printf("Query locus chr3:100: %d\n", query_locus(set, s1, 100));
    printf("Query locus chr4:100: %d\n", query_locus(set, s2, 100));
    printf("Query locus chr5:550: %d\n", query_locus(set, s3, 550));
    printf("Query overlap chr1:900-2200: %d\n", query_overlap(set, s4, 900, 2200));
    printf("Query overlap chr1:150-250: %d\n", query_overlap(set, s4, 150, 250));
    printf("Query overlap chr1:50-150: %d\n", query_overlap(set, s4, 50, 150));
    printf("Query overlap chr1:110-150: %d\n", query_overlap(set, s4, 110, 150));
    printf("Query overlap chr1:3000-4000: %d\n", query_overlap(set, s4, 3000, 4000));
    printf("Query overlap chr1:1-2: %d\n", query_overlap(set, s4, 1, 2));
    printf("\n");
    printf("Free Regions\n");
    free((char*) ks_release(s)); free(s);
    free((char*) ks_release(s1)); free(s1);
    free((char*) ks_release(s2)); free(s2);
    free((char*) ks_release(s3)); free(s3);
    free((char*) ks_release(s4)); free(s4);
    destroy_region_set(set);
    printf("Create Invalid Region:\n");
    s = (kstring_t*) calloc(1, sizeof(kstring_t));
    kputs("chr1:100-0assda0,", s);
    set = init_region_set(s, false);
    printf("set: %lu\n", (unsigned long) set);
    free((char*) ks_release(s)); free(s);
    destroy_region_set(set);
}
*/