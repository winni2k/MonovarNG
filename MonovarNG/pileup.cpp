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
#include "wrdouble.hpp"

#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <array>
#include <cmath>
#include <chrono>

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


void Pileup::setCombi(const Combination* combiPtr) {
    // sets combi
    combi = combiPtr;
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
    numCells = cells.size(); // set numcells to be the number of cells with read
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

vector<array<wrdouble, 3>> Pileup::computeLikelihoods(const array<array<array<double, 4>, 4>, 4>& genotypePriors, double pDropout) {
    // computes likelihoods L(g=0, 1, 2) for each cell
    vector<array<wrdouble, 3>> likelihoods;
    for (auto &cell: cells) {
        // g = 0
        wrdouble g0 = 1.0;
        for (int i = 0; i < cell.numReads; i++) {
            double probRead = genotypePriors[refBase][refBase][cell.bases[i]]; // likelihood of read given both refbase
            g0 *= cell.qualities[i]*(1-probRead)/3 + (1-cell.qualities[i])*probRead;
        }
        
        // g = 2
        wrdouble g2 = 1.0;
        for (int i = 0; i < cell.numReads; i++) {
            double probRead = genotypePriors[altBase][altBase][cell.bases[i]]; // likelihood of read given both altbase
            g2 *= cell.qualities[i]*(1-probRead)/3 + (1-cell.qualities[i])*probRead; 
        }
        
        // g = 1
        wrdouble g1 = 0, probADO = (g0+g2)/2.0, probNoADO = 1.0;
        for (int i = 0; i < cell.numReads; i++) {
            double probRead = genotypePriors[refBase][altBase][cell.bases[i]]; // likelihood of read given refbase & altbase
            probNoADO *= cell.qualities[i]*(1-probRead)/3 + (1-cell.qualities[i])*probRead; 
        }
        g1 = probADO * pDropout + probNoADO * (1-pDropout);
        
        likelihoods.push_back(array<wrdouble, 3>{g0, g1, g2});
    }
    
    return likelihoods;
}

vector<wrdouble> Pileup::computeDP(const vector<array<wrdouble, 3>>& likelihoods) {
    // computes dp for h_j,l and returns the row for j = numCells
    vector<vector<wrdouble>> mem;
    mem.resize(numCells); // j = [0, numCells) - 0-indexed
    for (int j = 0; j < numCells; j++) mem[j].resize(2*numCells+1, wrdouble(0)); // row j: l = [0, 2j+2]
    
    // Base case
    mem[0][0] = likelihoods[0][0];
    mem[0][1] = likelihoods[0][1] * 2.0;
    mem[0][2] = likelihoods[0][2];
    
    // Recursive cases
    for (int j = 1; j < numCells; j++) {
        mem[j][0] = mem[j-1][0]*likelihoods[j][0];
        mem[j][1] = mem[j-1][1]*likelihoods[j][0] + mem[j-1][0]*likelihoods[j][1]*2.0;
        for (int l = 2; l <= 2*j+2; l++) {
            mem[j][l] = mem[j-1][l]*likelihoods[j][0] + mem[j-1][l-1]*likelihoods[j][1]*2.0 + mem[j-1][l-2]*likelihoods[j][2];
        }
    }
    
    return mem[numCells-1]; // return the row for j = numCells
}

vector<wrdouble> Pileup::computeAltLikelihoods(const vector<wrdouble>& dp) {
    // computes alt count likelihoods, dividing each element i by 2*numCells C i
    vector<wrdouble> combis = combi->getRow(2*numCells);
    vector<wrdouble> altLikelihoods = dp;
    for (int i = 0; i <= 2*numCells; i++) altLikelihoods[i] /= combis[i];
    return altLikelihoods;
}

double Pileup::computeZeroVarProb(const array<array<array<double, 4>, 4>, 4>& genotypePriors, double pDropout) {
    // Generate variant number prior array
    vector<double> altCountPriors = genAltCountPriors(cellsWithRead());
    
    // Generate likelihoods L(g=0, 1, 2) for each cell
    vector<array<wrdouble, 3>> likelihoods = computeLikelihoods(genotypePriors, pDropout);
    
    // Generate dp
    vector<wrdouble> dp = computeDP(likelihoods);

    // Generate alternate count likelihoods
    vector<wrdouble> altLikelihoods = computeAltLikelihoods(dp);
    
    // Compute probability of mutation
    wrdouble base = 0;
    for (int i = 0; i <= 2*numCells; i++) {
        base += altLikelihoods[i] * altCountPriors[i];
    }
    wrdouble probability = altLikelihoods[0] * altCountPriors[0] / base; 
    
    return double(probability);
}

