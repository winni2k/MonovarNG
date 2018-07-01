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
    char altBase; // alternate base at position
    vector<SingleCellPos> cells; // data for individual cell reads
    vector<SingleCellPos> allCells; // data for all cells, an archived version of cells
    
    Pileup(int numCells, string row);
    
    void print(string filename = ""); // prints bases and qualities for debugging, and appends to file if specified
    
    int totalDepth(); // gets total depth (no. of reads)
    int refDepth(); // gets number of reads matching reference base
    
    int cellsWithRead(); // gets number of cells with reads
    int cellsWithAlt(); // gets number of cells with alternate alleles
    
    void sanitizeBases(); // Removes ins/deletions, special symbols, and cleans up all bases. Also changes refbase to upper
    void computeQualities(); // Converts the quality score string into decimal scores
    void filterCellsWithRead(); // archives cells to allCells, and filters cells for only those with reads
    
    array<int, 4> baseFreq(); // gets frequencies of each base - A, C, T, G
    bool setAltBase(); // sets the alternate base for position
    
    void convertBasesToInt(); // Converts all bases to integers: A=0, C=1, T=2, G=3, without changing data structure. Acts on cells and refbase/altbase

    
    

};

#endif /* pileup_hpp */
