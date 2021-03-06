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
#include "phred.hpp"
#include "ap.h"
#include "statistics.h"

#include <boost/algorithm/string.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <array>
#include <cmath>
#include <chrono>

using namespace std;
using namespace utility;

Pileup::Pileup(int numCells, string& row) : numCells(numCells) {
    // Parses row and saves them into tokens
    vector<string> tokens;
    tokens.reserve(3*(numCells+1));
    boost::split(tokens, row, boost::is_any_of("\t"));
    
    seqID = tokens[0];
    seqPos = stoi(tokens[1]);
    boost::trim(tokens[2]);
    refBase = toupper(tokens[2][0]);
    
    cells.reserve(numCells);
    for (int i = 0; i < numCells; i++) {
        cells.push_back(SingleCellPos(stoi(tokens[3*i+3]), tokens[3*i+4], tokens[3*i+5]));
    }
}

void Pileup::print(string filename, bool quality) {
    // prints bases and qualities for debugging
    if (filename.size()) freopen(filename.c_str(), "a", stdout); // save to file
    //    printf("%d\n", seqPos);
    printf("%d cells.\n", numCells);
    printf("Ref = %d, alt = %d\n", refBase, altBase);
    printf("Bases:\n");
    for (auto &cell: cells) {
        for (char c: cell.bases) printf("%d", c);
        printf("\n");
    }
    if (quality) {
        printf("Qualities:\n");
        for (auto &cell: cells) {
            for (double v: cell.qualities) printf("%lf\t", v);
            printf("\n");
        }
    }
}

void Pileup::setObjs(const Combination* combiPtr, const Phred* phredPtr) {
    // sets combi and phred
    combi = combiPtr;
    phred = phredPtr;
}


int Pileup::totalDepth() {
    int count = 0;
    for (SingleCellPos& cell: cells) {
        count += cell.numReads;
    }
    return count;
}

vector<pair<int, int>> Pileup::cellDepths() {
    // gets depth for each cell
    vector<pair<int, int>> depths;
    for (auto& cell: allCells) {
        int refCount = cell.countAllele(refBase), altCount = cell.countAllele(altBase);
        depths.push_back(make_pair(refCount, altCount));
    }
    return depths;
}

int Pileup::refDepth() {
    int count = 0;
    for (SingleCellPos& cell: cells) {
        count += cell.refCount();
    }
    return count;
}


int Pileup::cellsWithRead() {
    // gets number of cells with reads
    int count = 0;
    for (SingleCellPos& cell: cells) {
        count += cell.hasReads();
    }
    return count;
}

int Pileup::cellsWithAlt() {
    // gets number of cells with alternate alleles
    int count = 0;
    for (SingleCellPos& cell: cells) {
        count += cell.hasAltAllele();
    }
    return count;
}

void Pileup::filterCellsWithRead() { 
    // archives cells to allCells, and filters cells for only those with reads
    allCells = cells;
    cells = vector<SingleCellPos>();
    cells.reserve(allCells.size());
    for (auto& cell: allCells) {
        if (cell.hasReads()) cells.push_back(cell);
    }
    numCells = cells.size(); // set numcells to be the number of cells with read
}


void Pileup::sanitizeBases() {
    // Removes ins/deletions, special symbols, and cleans up all bases. Also changes refbase to upper. Also converts to numbers. Returns the number of forward and backward strands for each base.
    refBase = toupper(refBase);
    if (refBase == 'A') refBase = 0;
    else if (refBase == 'C') refBase = 1;
    else if (refBase == 'T') refBase = 2;
    else if (refBase == 'G') refBase = 3;
    
    for (int i = 0; i < 4; i++) for (int j = 0; j < 2; j++) strandCount[i][j] = 0;
    for (SingleCellPos &cell: cells) {
        auto strands = cell.sanitizeBases(refBase);
        for (int i = 0; i < 4; i++) for (int j = 0; j < 2; j++) strandCount[i][j] += strands[i][j];
        cell.truncateReads();
    }
}

void Pileup::computeQualities() {
    // Converts the quality score string into decimal scores
    for (auto &cell: cells) {
        cell.computeQuality(phred);
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
    int maxfreq = 0;
    for (int i = 0; i < 4; i++) {
        if (i != refBase && baseFrequency[i] >= maxfreq) {
            altBase = i;
            maxfreq = max(maxfreq, baseFrequency[i]);
        }
    }
    return maxfreq; // true if maxfreq > 0
}

vector<array<wrdouble, 3>> Pileup::computeLikelihoods(const array<array<array<double, 4>, 4>, 4>& genotypePriors, double pDropout) {
    // computes likelihoods L(g=0, 1, 2) for each cell
    vector<array<wrdouble, 3>> likelihoods;
    for (auto &cell: cells) {
        double collect0 = 1.0, collect1 = 1.0, collect2 = 1.0; // collects the products
        
        wrdouble g0 = 1.0, g2 = 1.0, probNoADO = 1.0;
        for (int i = 0; i < cell.numReads; i++) {
            // g = 0
            double probRead0 = genotypePriors[refBase][refBase][cell.bases[i]]; // likelihood of read given both refbase
            collect0 *= cell.qualities[i]*(1-probRead0)/3 + (1-cell.qualities[i])*probRead0;
            if (!(i%30)) {
                g0 *= collect0;
                collect0 = 1.0;
            }
            // g = 2
            double probRead2 = genotypePriors[altBase][altBase][cell.bases[i]]; // likelihood of read given both refbase
            collect2 *= cell.qualities[i]*(1-probRead2)/3 + (1-cell.qualities[i])*probRead2;
            if (!(i%30)) {
                g2 *= collect2;
                collect2 = 1.0;
            }
            // g = 1
            double probRead1 = genotypePriors[refBase][altBase][cell.bases[i]]; // likelihood of read given both refbase
            collect1 *= cell.qualities[i]*(1-probRead1)/3 + (1-cell.qualities[i])*probRead1;
            if (!(i%30)) {
                probNoADO *= collect1;
                collect1 = 1.0;
            }
        }
        g0 *= collect0;
        g2 *= collect2;
        probNoADO *= collect1;
        
        wrdouble probADO = (g0+g2)/2.0;
        wrdouble g1 = probADO * pDropout + probNoADO * (1-pDropout);
        
        likelihoods.push_back(array<wrdouble, 3>{g0, g1, g2});
    }
    
    return likelihoods;
}

vector<wrdouble> Pileup::computeDP(const vector<array<wrdouble, 3>>& likelihoods) {
    // computes dp for h_j,l and returns the row for j = numCells
    vector<vector<wrdouble>> mem;
    mem.resize(likelihoods.size()); // j = [0, numCells) - 0-indexed
    for (int j = 0; j < likelihoods.size(); j++) mem[j].resize(2*likelihoods.size()+1, wrdouble(0)); // row j: l = [0, 2j+2]
    wrdouble wr2 = 2.0;
    
    // Base case
    mem[0][0] = likelihoods[0][0];
    mem[0][1] = likelihoods[0][1] * wr2;
    mem[0][2] = likelihoods[0][2];
    
    // Recursive cases
    for (int j = 1; j < likelihoods.size(); j++) {
        mem[j][0] = mem[j-1][0]*likelihoods[j][0];
        mem[j][1] = mem[j-1][1]*likelihoods[j][0] + mem[j-1][0]*likelihoods[j][1]*wr2;
        for (int l = 2; l <= 2*j+2; l++) {
            mem[j][l] = mem[j-1][l]*likelihoods[j][0] + mem[j-1][l-1]*likelihoods[j][1]*wr2 + mem[j-1][l-2]*likelihoods[j][2];
        }
    }
    
    return mem[likelihoods.size()-1]; // return the row for j = numCells
}

vector<wrdouble> Pileup::computeAltLikelihoods(const vector<wrdouble>& dp) {
    // computes alt count likelihoods, dividing each element i by 2*numCells C i
    vector<wrdouble> combis = combi->getRow(2*numCells);
    vector<wrdouble> altLikelihoods = dp;
    for (int i = 0; i <= 2*numCells; i++) altLikelihoods[i] /= combis[i];
    return altLikelihoods;
}

wrdouble Pileup::computeZeroVarProb(const array<array<array<double, 4>, 4>, 4>& genotypePriors, double pDropout) {
    // Generate variant number prior array
    vector<double> altCountPriors = genAltCountPriors(cellsWithRead());
    
    // Generate likelihoods L(g=0, 1, 2) for each cell
    likelihoodsGlob = computeLikelihoods(genotypePriors, pDropout);
//    printf("Likelihoods:\n");
//    for (int i = 0; i < numCells; i++) {
//        for (int j = 0; j < 3; j++) cout << likelihoodsGlob[i][j] << "\t";
//        printf("\n");
//    }
    // Generate dp
    vector<wrdouble> dp = computeDP(likelihoodsGlob);
//    printf("DP:\n");
//    for (int i = 0; i < numCells*2+1; i++) cout << dp[i] << "\t";
//    cout << endl;
    
    // Generate alternate count likelihoods
    vector<wrdouble> altLikelihoods = computeAltLikelihoods(dp);
//    printf("Alt likelihoods:\n");
//    for (int i = 0; i < numCells*2+1; i++) cout << altLikelihoods[i] << "\t";
//    cout << endl;
    
    // Compute probability of mutation
    probBase = 0.0;
    for (int i = 0; i <= 2*numCells; i++) {
        probBase += altLikelihoods[i] * altCountPriors[i];
    }
    wrdouble probability = (altLikelihoods[0] * altCountPriors[0]) / probBase; 
    return probability;
}

double Pileup::computeC(int l, int v) {
    // computes the C function
    return combi->getValue(l, v) * combi->getValue(2*numCells-l, 2-v) / combi->getValue(2*numCells, 2);
//    if (v == 0) {
//        return (2.0*numCells-l)*(2.0*numCells-l-1)/((2.0*numCells)*(2.0*numCells-1));
//    } else if (v == 1) {
//        return 2.0*l*(2*numCells-l)/((2.0*numCells)*(2.0*numCells-1));
//    } else {
//        return (1.0*l*(l-1))/((2.0*numCells)*(2.0*numCells-1));
//    }
}

vector<int> Pileup::computeGenotype() {
    // computes the genotype of each cell, 0, 1 or 2
    vector<double> altCountPriors = genAltCountPriors(numCells);
//    printf("Alt count priors:\n");
//    for (int i = 0; i <= numCells*2; i++) printf("%lf\t", altCountPriors[i]);
//    printf("\n");
    
//    cout << "ProbBase = " << probBase << endl;
    vector<int> genotypes; genotypes.reserve(numCells);
    for (int i = 0; i < numCells; i++) {
        // Build new likelihoods, by removing this cell
        vector<array<wrdouble, 3>> newLikelihoods;
        newLikelihoods.reserve(numCells-1);
        for (int j = 0; j < numCells; j++) {
            if (j != i) newLikelihoods.push_back(likelihoodsGlob[j]);
        }
        vector<wrdouble> dp;
        if (numCells != 1) dp = computeDP(newLikelihoods);
        
//        printf("\nCell %d\n", i);
//        printf("DP:\n");
//        for (wrdouble j: dp) cout << j << "\t";
//        cout << endl;
        
        wrdouble probs[3]; // probability of each genotype
        for (int j = 0; j < 3; j++) {
            probs[j] = 0.0;
            if (numCells != 1) {
                for (int l = j; l <= 2*numCells-2+j; l++) {
    //                cout << "Base: " << probBase << endl;
                    probs[j] += dp[l-j]*((computeC(l, j)*altCountPriors[l]));
//                    cout << dp[l-j] << " " << computeC( l, j) << " " << altCountPriors[l] << endl; 
//                    cout << "prob: " << probs[j] << endl;
                }
            } else probs[j] = altCountPriors[i]; // There aren't any other cells
            
//            cout << "before multiply" << probs[j] << "\t";
            probs[j] *= likelihoodsGlob[i][j];
//            cout << "after multiply" << probs[j] << "\t";
            probs[j] /= probBase;
//            cout << probs[j] << "\t";
//            cout << "finalprob: " << probs[j] << "\n";
        }
//        cout << endl;
        
        int bestGenotype = -1;
        wrdouble highestProb = 0;
        for (int j = 0; j < 3; j++) {
            if (probs[j] > highestProb) {
                highestProb = probs[j];
                bestGenotype = j;
            }
        }
        genotypes.push_back(bestGenotype);
    }
    
    // Add '-1' for cells with no reads
    vector<int> allGenotypes;
    int pos = 0; // position in genotypes
    for (auto& cell: allCells) {
        if (cell.hasReads()) {
            allGenotypes.push_back(genotypes[pos]);
            pos++;
        } else allGenotypes.push_back(-1);
    }
    
    if (pos != genotypes.size()) exit(1); // error!
    
    return allGenotypes;
}


double Pileup::computeWilcoxon() {
    // computes the Mann-Whitney-Wilcoxon test
    vector<double> ref, alt;
    for (auto& cell: cells) {
        for (int i = 0; i < cell.numReads; i++) {
            if (cell.bases[i] == refBase) ref.push_back(cell.qualities[i]);
            else if (cell.bases[i] == altBase) alt.push_back(cell.qualities[i]);
        }
    }
    if (ref.size() < 5 || alt.size() < 5) return 0.0;
    alglib::real_1d_array refA, altA;
    refA.setcontent(ref.size(), ref.data());
    altA.setcontent(alt.size(), alt.data());
    double bothtails, lefttail, righttail;
    alglib::mannwhitneyutest(altA, alt.size(), refA, ref.size(), bothtails, lefttail, righttail);
    return bothtails;
}

double Pileup::qualityByDepth(const double& quality, const vector<int>& genotype) {
    // computes QualByDepth, quality divided by the number of reads in cells with mutation
    int depth = 0;
    for (int i = 0; i < allCells.size(); i++) {
        if (genotype[i] == 1 || genotype[i] == 2) {
            depth += allCells[i].numReads;
        }
    }
    if (depth == 0) return quality;
    else return quality/depth;
}

double Pileup::computeStrandBias() {
    // computes strand bias
    for (int i = 0; i < 4; i++) for (int j = 0; j < 2; j++) strandCount[i][j]++;
    
    double r = double(strandCount[0][0]*strandCount[1][1])/(strandCount[0][1]*strandCount[1][0]);
    double refRatio = double(min(strandCount[0][0], strandCount[0][1]))/max(strandCount[0][0], strandCount[0][1]);
    double altRatio = double(min(strandCount[1][0], strandCount[1][1]))/max(strandCount[1][0], strandCount[1][1]);
    
    return log(refRatio/altRatio*(r+1.0/r));
}

double Pileup::psarr(vector<pair<int, int>>& depths) {
    // computes PSARR, ratio of per-sample alt allele to ref allele
    int refCount = 0, altCount = 0, refReads = 0, altReads = 0;
    for (auto& counts: depths) {
        if (counts.first + counts.second == 0) continue;
        if (counts.first) refReads++;
        if (counts.second) altReads++;
        refCount += counts.first;
        altCount += counts.second;
    }
    if (altReads == 0) altReads = 1.0;
    if (refReads == 0) refReads = 1.0;
    double num = double(altCount)/altReads;
    double den = double(refCount)/refReads;
    if (den == 0) den = 1.0;
    return num/den;
}








