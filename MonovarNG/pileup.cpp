//
//  pileup.cpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "pileup.hpp"
#include "single_cell_pos.hpp"

#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <iostream>

using namespace std;

Pileup::Pileup(int numCells, string row) : numCells(numCells) {
    boost::trim(row);
    vector<string> tokens;
    boost::split(tokens, row, boost::is_any_of("\t"));
    
    seqID = tokens[0];
    seqPos = stoi(tokens[1]);
    boost::trim(tokens[2]);
    refBase = toupper(tokens[2][0]);
    for (int i = 0; i < numCells; i++) {
        cells.push_back(SingleCellPos(stoi(tokens[3*i+3]), tokens[3*i+4], tokens[3*i+5], refBase));
    }
}

int Pileup::totalDepth() {
    int count = 0;
    for (SingleCellPos cell: cells) {
        count += cell.numReads;
    }
    return count;
}

int Pileup::refDepth() {
    int count = 0;
    for (SingleCellPos cell: cells) {
        count += cell.refCount();
    }
    return count;
}

int Pileup::cellsWithRead() {
    // gets number of cells with reads
    int count = 0;
    for (SingleCellPos cell: cells) {
        count += cell.hasReads();
    }
    return count;
}
int Pileup::cellsWithAlt() {
    // gets number of cells with alternate alleles
    int count = 0;
    for (SingleCellPos cell: cells) {
        count += cell.hasAltAllele();
    }
    return count;
}









