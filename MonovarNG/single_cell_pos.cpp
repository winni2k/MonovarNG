//
//  single_cell_pos.cpp
//  MonovarNG
//
//  Created by tianyi on 6/26/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "single_cell_pos.hpp"

#include <string>
#include <vector>

using namespace std;

SingleCellPos::SingleCellPos(int numReads, string bases, string qualities, char refBase) : numReads(numReads), bases(bases), qualities(qualities), refBase(refBase) {} // initializer sets default values

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
