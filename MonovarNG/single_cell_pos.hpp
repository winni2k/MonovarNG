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
    string qualities; // quality at position, for cell
    char refBase; // reference base at position
    
    SingleCellPos(int numReads, string bases, string qualities, char refBase);
    
    int refCount(); // gets number of forward + backward matching reads matching reference base
    bool hasReads(); // gets whether the cell has reads (nonzero read depth)
    bool hasAltAllele(); // gets whether the cell has alternate alleles (non ,.)
};

#endif /* single_cell_pos_hpp */
