//
//  pileup.cpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "utility.hpp"
#include "pileup.hpp"
#include "single_cell_pos.hpp"

#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <array>
#include <cmath>

using namespace std;
using namespace utility;

Pileup::Pileup(int numCells, string row) : numCells(numCells) {
    boost::trim(row);
    vector<string> tokens;
    boost::split(tokens, row, boost::is_any_of("\t"));
    
    seqID = tokens[0];
    seqPos = stoi(tokens[1]);
    boost::trim(tokens[2]);
    refBase = toupper(tokens[2][0]);
    for (int i = 0; i < numCells; i++) {
        cells.push_back(SingleCellPos(stoi(tokens[3*i+3]), tokens[3*i+4], tokens[3*i+5]));
    }
}

void Pileup::print(string filename) {
    // prints bases and qualities for debugging
    if (filename.size()) freopen(filename.c_str(), "a", stdout); // save to file
    //    printf("%d\n", seqPos);
    for (auto &cell: cells) {
        for (char c: cell.bases) printf("%d", c);
        printf("\n");
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

void Pileup::filterCellsWithRead() { 
    // archives cells to allCells, and filters cells for only those with reads
    allCells = cells;
    cells = vector<SingleCellPos>();
    for (auto cell: allCells) {
        if (cell.hasReads()) cells.push_back(cell);
    }
}


void Pileup::sanitizeBases() {
    // Removes ins/deletions, special symbols, and cleans up all bases. Also changes refbase to upper
    refBase = toupper(refBase);
    for (SingleCellPos &cell: cells) {
        cell.removeInsDels();
        cell.removeStartEnd();
        cell.cleanupBases(refBase);
        cell.truncateReads();
    }
}

void Pileup::computeQualities() {
    // Converts the quality score string into decimal scores
    for (auto &cell: cells) {
        cell.computeQuality();
    }
}

array<int, 4> Pileup::baseFreq() {
    // gets frequencies of each base - A, C, T, G
    array<int, 4> freq = {0};
    for (auto cell: cells) {
        array<int, 4> cellFreq = cell.baseFreq();
        for (int i = 0; i < 4; i++) freq[i] += cellFreq[i];
    }
    return freq;
}

bool Pileup::setAltBase() {
    // sets the alternate base for position, returning true for successful set
    array<int, 4> baseFrequency = baseFreq(); // A, C, T, G
    char corrBases [4] = {'A', 'C', 'T', 'G'}; // bases corresponding to position
    int maxfreq = 0;
    for (int i = 0; i < 4; i++) {
        if (corrBases[i] != refBase && baseFrequency[i] >= maxfreq) {
            altBase = corrBases[i];
            maxfreq = max(maxfreq, baseFrequency[i]);
        }
    }
    return maxfreq; // true if maxfreq > 0
}

void Pileup::convertBasesToInt() {
    // Converts all bases to integers: A=0, C=1, T=2, G=3, without changing data structure. Acts on cells and refbase/altbase
    string mapping = "ACTG";
    for (int i = 0; i < 4; i++) {
        if (refBase == mapping[i]) refBase = i;
        if (altBase == mapping[i]) altBase = i;
    }
    
    for (auto &cell: cells) cell.convertBasesToInt();
}

vector<array<double, 3>> Pileup::computeLikelihoods(array<array<array<double, 4>, 4>, 4> genotypePriors, double pDropout) {
    // computes likelihoods L(g=0, 1, 2) for each cell
    vector<array<double, 3>> likelihoods;
    for (auto &cell: cells) {
        // g = 0
        double g0 = 1.0;
        for (int i = 0; i < cell.numReads; i++) {
            double probRead = genotypePriors[refBase][refBase][cell.bases[i]]; // likelihood of read given both refbase
            g0 *= cell.qualities[i]*(1-probRead)/3 + (1-cell.qualities[i])*probRead; 
        }
        
        // g = 2
        double g2 = 1.0;
        for (int i = 0; i < cell.numReads; i++) {
            double probRead = genotypePriors[altBase][altBase][cell.bases[i]]; // likelihood of read given both altbase
            g2 *= cell.qualities[i]*(1-probRead)/3 + (1-cell.qualities[i])*probRead; 
        }
        
        // g = 1
        double g1 = 0, probADO = 1.0, probNoADO = (g0+g2)/2;
        for (int i = 0; i < cell.numReads; i++) {
            double probRead = genotypePriors[refBase][altBase][cell.bases[i]]; // likelihood of read given refbase & altbase
            probADO *= cell.qualities[i]*(1-probRead)/3 + (1-cell.qualities[i])*probRead; 
        }
        g1 = pDropout * probADO + (1-pDropout) * probNoADO;
        
        printf("%d %d ", refBase, altBase);
        for (char c: cell.bases) {
            printf("%d", c);
        }
        printf(" %lf %lf %lf\n", log10(g0), log10(g1), log10(g2));
        likelihoods.push_back(array<double, 3>{g0, g1, g2});
    }
    
    return likelihoods;
}


double Pileup::computeZeroVarProb(array<array<array<double, 4>, 4>, 4> genotypePriors, double pDropout) {
    // Generate variant number prior array
    vector<double> altCountPriors = genAltCountPriors(cellsWithRead());
    // Generate likelihoods L(g=0, 1, 2) for each cell
    vector<array<double, 3>> likelihoods = computeLikelihoods(genotypePriors, pDropout);
    
    return 0.0;
}

