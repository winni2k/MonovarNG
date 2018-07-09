//
//  single_cell_pos.cpp
//  MonovarNG
//
//  Created by tianyi on 6/26/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "single_cell_pos.hpp"
#include "phred.hpp"

#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

SingleCellPos::SingleCellPos(int numReads, string& bases, string& qualityString) : numReads(numReads), bases(bases), qualityString(qualityString) {} // initializer sets default values

int SingleCellPos::refCount() {
    // Returns number of forward + backward matching reads matching reference base
    int count = 0;
    for (char& base: bases) {
        if (base == '.' || base == ',') count++;
    }   
    return count;
}

bool SingleCellPos::hasReads() {
    // Returns whether the cell has read support (nonzero reads)
    return numReads;
}

bool SingleCellPos::hasAltAllele() {
    // Returns whether the cell has alternate alleles (non ,.)
    return numReads - refCount();
}

void SingleCellPos::sanitizeBases(char refBase) { 
    // remove ins/deletions, special symbols, and cleans up all bases. Also converts to numbers
    string newBases = "";
    newBases.reserve(bases.size());
    
    int state = 0; // 0 = normal, 1 = counting no. of ins/del, 2 = deleting ins/dels
    int baseCount = 0; // count of the number of bases inserted/deleted
    for (int i = 0; i < bases.size(); i++) {
        char& c = bases[i];
        if (state == 0) {
            if (c == '+' || c == '-') {
                // Insertion/deletion begins
                state = 1;
                baseCount = 0;
            } else {
                // Run other filters - start/end, and sanitize
                if (c == '^') { // start/end
                    i++; // skip one more character
                } else if (c != '$') {
                    // sanitization
                    if (c == '.' || c == ',' || c == '*') newBases += refBase;
                    else if (c == 'A' || c == 'a') newBases += char(0);
                    else if (c == 'C' || c == 'c') newBases += char(1);
                    else if (c == 'T' || c == 't') newBases += char(2);
                    else if (c == 'G' || c == 'g') newBases += char(3);
                }
            }
        } else if (state == 1) {
            if (isdigit(c)) {
                baseCount *= 10;
                baseCount += (c - '0');
            } else if (baseCount){
                // Is base, so start removing
                // Test for whether basecount is nonzero, so that we don't 
                state = 2;
                baseCount--;
                if (baseCount == 0) state = 0;
            } else {
                state = 0;
                // Run other filters - start/end, and sanitize
                if (c == '^') { // start/end
                    i++; // skip one more character
                } else if (c != '$') {
                    // sanitization
                    if (c == '.' || c == ',' || c == '*') newBases += refBase;
                    else if (c == 'A' || c == 'a') newBases += char(0);
                    else if (c == 'C' || c == 'c') newBases += char(1);
                    else if (c == 'T' || c == 't') newBases += char(2);
                    else if (c == 'G' || c == 'g') newBases += char(3);
                }
            }
        } else if (state == 2) {
            baseCount--;
            if (baseCount == 0) state = 0;
        }
    }
    
    bases = newBases;
}

void SingleCellPos::truncateReads() {
    // truncates numReads, bases and qualities to the shortest length; a naive way of dealing with input deviations
    int minLength = min({numReads, (int)bases.size(), (int)qualityString.size()});
    numReads = minLength;
    bases.resize(minLength);
    qualityString.resize(minLength);
}

void SingleCellPos::computeQuality(const Phred* phred) {
    qualities.reserve(numReads);
    // Converts the quality score string into decimal scores
    for (char& c: qualityString) {
        int phredVal = (int) c - 33;
        qualities.push_back(phred->qualities[phredVal]);
    }
}


array<int, 4> SingleCellPos::baseFreq() {
    // gets frequencies of each base - A, C, T, G
    array<int, 4> freq = {0};
    for (char& c: bases) {
        freq[c]++;
    }
    return freq;
}

