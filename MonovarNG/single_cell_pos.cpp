//
//  single_cell_pos.cpp
//  MonovarNG
//
//  Created by tianyi on 6/26/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "single_cell_pos.hpp"

#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

SingleCellPos::SingleCellPos(int numReads, string bases, string qualityString) : numReads(numReads), bases(bases), qualityString(qualityString) {} // initializer sets default values

int SingleCellPos::refCount() {
    // Returns number of forward + backward matching reads matching reference base
    int count = 0;
    for (char base: bases) {
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

void SingleCellPos::removeInsDels() { 
    // removes all insertions and deletions from the bases
    string newBases = "";
    int state = 0; // 0 = normal, 1 = counting no. of ins/del, 2 = deleting ins/dels
    int baseCount = 0; // count of the number of bases inserted/deleted
    for (char c: bases) {
        if (state == 0) {
            if (c == '+' || c == '-') {
                // Insertion/deletion begins
                state = 1;
                baseCount = 0;
            } else {
                newBases += c;
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
                newBases += c;
            }
        } else if (state == 2) {
            baseCount--;
            if (baseCount == 0) state = 0;
        }
    }
    
    bases = newBases;
}

void SingleCellPos::removeStartEnd() {
    // removes all start and end read symbols
    // note that a ^ is followed by another character denoting quality
    string newBases = "";
    for (int i = 0; i < bases.size(); i++) {
        char c = bases[i];
        if (c != '$' && c != '^') {
            newBases += c;
        } else if (c == '^') {
            i++; // skip one more character
        }
    }
    bases = newBases;
}

void SingleCellPos::cleanupBases(char refBase) {
    // cleans up bases such that it contains 'ACTG' only
    string newBases = "";
    for (char c: bases) {
        if (string(",.*actgACTG").find(c) != string::npos) {
            if (c == '.' || c == ',' || c == '*') newBases += refBase;
            else newBases += toupper(c);
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

void SingleCellPos::computeQuality() {
    // Converts the quality score string into decimal scores
    for (char c: qualityString) {
        int phred = (int) c - 33;
        qualities.push_back(pow((double) 10.0, -(double)phred/10));
    }
}

void SingleCellPos::convertBasesToInt() {
    // Converts all bases to integers: A=0, C=1, T=2, G=3, without changing data structure
    for (int i = 0; i < bases.size(); i++) {
        string mapping = "ACTG";
        for (int j = 0; j < 4; j++) {
            if (bases[i] == mapping[j]) bases[i] = j;
        }
    }
}


array<int, 4> SingleCellPos::baseFreq() {
    // gets frequencies of each base - A, C, T, G
    array<int, 4> freq = {0};
    for (char c: bases) {
        switch (c) {
            case 'A':
                freq[0]++;
                break;
            case 'C':
                freq[1]++;
                break;
            case 'T':
                freq[2]++;
                break;
            case 'G':
                freq[3]++;
                break;
        }
    }
    return freq;
}

