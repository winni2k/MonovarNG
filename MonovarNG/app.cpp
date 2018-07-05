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

#include <cstdio>
#include <iostream>
#include <array>

using namespace std;
using namespace utility;

App::App(Config config, vector<string> bamIDs, vector<Pileup> pileupRows) : mutationThreshold(config.mutationThreshold), pFalsePositive(config.pFalsePositive), pDropout(config.pDropout), numThreads(config.numThreads), useConsensusFilter(config.useConsensusFilter), positions(pileupRows), combi(Combination(2*bamIDs.size())) {
    numCells = bamIDs.size();
    
    // Write some VCF stuff
    output = VCFDocument();
    output.write_stuff(bamIDs);
    
    numPos = positions.size();
    
    // Set combi object for each pileup
    for (auto& row: positions) {
        row.setCombi(&combi);
    }
}

void App::runAlgo() {
    for (int rowN = 0; rowN < numPos; rowN++) {
        Pileup& position = positions[rowN];
                
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
            if ((rowN+1) % 1000 == 0) printf("Processed row %d\n", rowN+1);
            continue;
        }
        
        // Parse and filter reads
        position.sanitizeBases();
        position.filterCellsWithRead();
        position.computeQualities();
        
        // Another filtration, in case all reads are erased
        if (!position.cellsWithRead()) continue;
        
        // Find and set alternate base at positon. If alt base cannot be set, return
        if (!position.setAltBase()) continue;
        
        // Convert all bases to integers to speedup computation
        position.convertBasesToInt(); 
        
        // Generate genotype priors
        array<array<array<double, 4>, 4>, 4> genotypePriors; // probability of read given genotype e.g. P(^AA)(_C)
        if (position.cellsWithRead() > (numCells/2)-1 && position.cellsWithAlt() == 1) genotypePriors = genGenotypePriors(0.2);
        else if (position.cellsWithRead() > (numCells/2) && position.cellsWithAlt() == 2 && position.totalDepth() > 30 && altFreq < 0.1) genotypePriors = genGenotypePriors(0.1);
        else genotypePriors = genGenotypePriors(pFalsePositive);
        
        // Compute probability of zero mutations given data
        double zeroVarProb = position.computeZeroVarProb(genotypePriors, pDropout);
        if (zeroVarProb < 0.5) {
            printf("%d %lf\n", rowN+1, zeroVarProb);
        }
        
        
        
        
        
        
        
        if ((rowN+1) % 1000 == 0) printf("Processed row %d\n", rowN+1);
    }
}
