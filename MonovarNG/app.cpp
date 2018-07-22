//
//  app.cpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "app.hpp"
#include "vcf.hpp"
#include "utility.hpp"
#include "combination.hpp"
#include "wrdouble.hpp"
#include "ThreadPool.h"

#include <cstdio>
#include <iostream>
#include <array>
#include <mutex>

using namespace std;
using namespace utility;

App::App(Config& config, vector<string>& bamIDs, vector<string>& pileupRows) : mutationThreshold(config.mutationThreshold), pFalsePositive(config.pFalsePositive), pDropout(config.pDropout), numThreads(config.numThreads), useConsensusFilter(config.useConsensusFilter), pileup(pileupRows), combi(Combination(2*bamIDs.size())), phred(Phred()), output(VCFDocument(config.outputFilename)) {
    numCells = bamIDs.size();
    
    // Write some VCF stuff
    output.writeDefHeader();
    vector<string> bamFilenames = getBamFilenames(config.bamfileNames);
    output.writeHeaderInfo(config.referenceFilename, bamFilenames);
    
    numPos = pileup.size();
    
    // Set combi object for each pileup
//    for (auto& row: positions) {
//        row.setObjs(&combi, &phred);
//    }
}

void App::processRow(int rowN) {
    // processes row of data
    Pileup position = getPileup(numCells, pileup[rowN]);
    position.setObjs(&combi, &phred);
        
    int totalDepth = position.totalDepth(), refDepth = position.refDepth(); // total no. of reads / no. matching reference base
    
    int altCount = totalDepth - refDepth; // no. of alternate reads
    
    double altFreq;
    if (totalDepth == 0) altFreq = 0;
    else altFreq = (double) altCount / totalDepth;
    
    // Prefilter
    int prefilter = 0;
    if (totalDepth == refDepth) prefilter = 1; // no reads supporting alternate allele, so no operations needed
    else if (totalDepth > 30 && (altCount <= 2 || altFreq <= 0.001)) prefilter = 2; // prefiltered due to unlikely mutation
    else if (string("ATGC").find(position.refBase) == string::npos) prefilter = 3; // bad reference
    else if (totalDepth <= 10) prefilter = 4; // insufficient data
    if (prefilter) {
        //            if (prefilter == 3) printf("%d Prefiltered due to %d\n", rowN+1, prefilter);
        if ((rowN+1) % 50000 == 0) printf("Processed row %d\n", rowN+1);
        return;
    }
    
    // Parse and filter reads
    position.sanitizeBases();
    position.filterCellsWithRead();
    
    // Another filtration, in case all reads are erased
    if (!position.numCells) return;
    
    // Find and set alternate base at positon. If alt base cannot be set, return
    if (!position.setAltBase()) return;
    
    // Compute quality scores
    position.computeQualities();
    
    // Generate genotype priors
    array<array<array<double, 4>, 4>, 4> genotypePriors; // probability of read given genotype e.g. P(^AA)(_C)
    int cellsWithAlt = position.cellsWithAlt();
    if (position.cellsWithRead() > (numCells/2)-1 && cellsWithAlt == 1) genotypePriors = genGenotypePriors(0.2);
    else if (position.cellsWithRead() > (numCells/2) && cellsWithAlt == 2 && position.totalDepth() > 30 && altFreq < 0.1) genotypePriors = genGenotypePriors(0.1);
    else genotypePriors = genGenotypePriors(pFalsePositive);
    
    // Compute probability of zero mutations given data
    wrdouble zeroVarProb = position.computeZeroVarProb(genotypePriors, pDropout);
    if (zeroVarProb < 0.05) {
        //            cout << position.seqID << " " << position.seqPos << " " << zeroVarProb << endl;
        //            position.print("", true);
        vector<int> genotypes = position.computeGenotype();
        //            for (int i = 0; i < genotypes.size(); i++) cout << genotypes[i] << "\t";
        //            cout << endl;
        double quality = zeroVarProb.phred();
        double qualityByDepth = position.qualityByDepth(quality, genotypes);
        double strandBias = position.computeStrandBias();
        vector<pair<int, int>> cellDepths = position.cellDepths();
        double psarr = position.psarr(cellDepths);
        
        outputMutex.lock();
        output.writeRow(position.seqID, position.seqPos, position.refBase, position.altBase, quality, position.computeWilcoxon(), qualityByDepth, strandBias, psarr, genotypes, position.totalDepth(), cellDepths, position.likelihoodsGlob);
        outputMutex.unlock();
    }
    
    if ((rowN+1) % 50000 == 0) printf("Processed row %d\n", rowN+1);
}

void App::runAlgo() {
    ThreadPool pool(numThreads);
    for (int rowN = 0; rowN < numPos; rowN++) {
        pool.enqueue(&App::processRow, this, rowN);
//        processRow(rowN); // single threaded
    }
}
