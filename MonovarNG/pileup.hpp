//
//  pileup.hpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef pileup_hpp
#define pileup_hpp

#include "single_cell_pos.hpp"

#include <stdio.h>
#include <string>
#include <vector>

using namespace std;

struct Pileup {
    // Stores data in a row in pileup format
    int numCells;
    string seqID; // sequence identifier
    int seqPos; // position in sequence (starting from 1)
    char refBase; // reference base at position
    vector<SingleCellPos> cells; // data for individual cell reads 
    
    Pileup(int numCells, string row);
    
    int totalDepth(); // gets total depth (no. of reads)
    int refDepth(); // gets number of reads matching reference base
    
    int cellsWithRead(); // gets number of cells with reads
    int cellsWithAlt(); // gets number of cells with alternate alleles
};

#endif /* pileup_hpp */
