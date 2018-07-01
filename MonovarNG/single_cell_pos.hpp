//
//  single_cell_pos.hpp
//  MonovarNG
//
//  Created by tianyi on 6/26/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef single_cell_pos_hpp
#define single_cell_pos_hpp

#include <stdio.h>
#include <string>
#include <vector>

using namespace std;


struct SingleCellPos {
    int numReads; // number of reads at position, for cell
    string bases; // bases at position, for cell
    string qualityString; // quality string for qualities at position, for cell
    vector<double> qualities; // qualities for each read at position, for cell
    
    SingleCellPos(int numReads, string bases, string qualityString);
    
    int refCount(); // gets number of forward + backward matching reads matching reference base
    bool hasReads(); // gets whether the cell has reads (nonzero read depth)
    bool hasAltAllele(); // gets whether the cell has alternate alleles (non ,.)
    
    void removeInsDels(); // removes all insertions and deletions from the bases
    void removeStartEnd(); // removes all start and end read symbols
    void cleanupBases(char refBase); // cleans up bases such that it contains 'ACTG' only
    void truncateReads(); // truncates numReads, bases and qualities to the shortest length; a naive way of dealing with input deviations
    void computeQuality(); // Converts the quality score string into decimal scores
    void convertBasesToInt(); // Converts all bases to integers: A=0, C=1, T=2, G=3, without changing data structure
    
    array<int, 4> baseFreq(); // gets frequencies of each base - A, C, T, G
    
};

#endif /* single_cell_pos_hpp */
