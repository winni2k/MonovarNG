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
#include "wrdouble.hpp"
#include "combination.hpp"
#include "phred.hpp"

#include <stdio.h>
#include <string>
#include <vector>
#include <array>

using namespace std;

//class Combination; // forward definition of combination object

struct Pileup {
    // Stores data in a row in pileup format
    int numCells;
    string seqID; // sequence identifier
    int seqPos; // position in sequence (starting from 1)
    char refBase; // reference base at position
    char altBase; // alternate base at position
    vector<SingleCellPos> cells; // data for individual cell reads
    vector<SingleCellPos> allCells; // data for all cells, an archived version of cells
    
    const Combination* combi; // computes nCr, as a row of nC0...nCn
    const Phred* phred; // computes phred quality scores
    
    vector<array<wrdouble, 3>> likelihoodsGlob; // Likelihoods, saved from zeroVarProb for use in genotyping
    wrdouble probBase; // base, sum0_2m p(D|l)p(l) 
    
    Pileup(int numCells, string& row);
    
    void print(string filename = ""); // prints bases and qualities for debugging, and appends to file if specified
    
    void setObjs(const Combination* combiPtr, const Phred* phred); // sets combi and phred
    
    int totalDepth(); // gets total depth (no. of reads)
    int refDepth(); // gets number of reads matching reference base
    
    int cellsWithRead(); // gets number of cells with reads
    int cellsWithAlt(); // gets number of cells with alternate alleles
    
    void sanitizeBases(); // removes ins/deletions, special symbols, and cleans up all bases. Also changes refbase to upper
    void computeQualities(); // converts the quality score string into decimal scores
    void filterCellsWithRead(); // archives cells to allCells, and filters cells for only those with reads
    
    array<int, 4> baseFreq(); // gets frequencies of each base - A, C, T, G
    bool setAltBase(); // sets the alternate base for position
    
    void convertBasesToInt(); // converts all bases to integers: A=0, C=1, T=2, G=3, without changing data structure. Acts on cells and refbase/altbase

    
    vector<array<wrdouble, 3>> computeLikelihoods(const array<array<array<double, 4>, 4>, 4>& genotypePriors, double pDropout); // computes likelihoods L(g=0, 1, 2) for each cell
    vector<wrdouble> computeDP(const vector<array<wrdouble, 3>>& likelihoods); // computes dp for h_j,l and returns the row for j = numCells
    vector<wrdouble> computeAltLikelihoods(const vector<wrdouble>& dp); // computes alt count likelihoods, dividing each element i by 2*numCells C i
    wrdouble computeZeroVarProb(const array<array<array<double, 4>, 4>, 4>& genotypePriors, double pDropout); // computes the probability of zero mutations given data
    double computeC(int l, int v); // computes the C function
    vector<int> computeGenotype(); // computes the genotype of each cell, 0, 1 or 2
};

#endif /* pileup_hpp */
