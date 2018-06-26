//
//  app.cpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "app.hpp"
#include "vcf.hpp"

#include <cstdio>
#include <iostream>

App::App(Config config, vector<string> bamIDs, vector<Pileup> pileupRows) : mutationThreshold(config.mutationThreshold), pFalsePositive(config.pFalsePositive), pDropout(config.pDropout), numThreads(config.numThreads), useConsensusFilter(config.useConsensusFilter), positions(pileupRows) {
    
    numCells = bamIDs.size();
    
    // Write some VCF stuff
    output = VCFDocument();
    output.write_stuff(bamIDs);
    
    numPos = positions.size();
}

void App::runAlgo() {
    for (int rowN = 0; rowN < numPos; rowN++) {
        Pileup position = positions[rowN];
        
        int totalDepth = position.totalDepth(), refDepth = position.refDepth(); // total no. of reads / no. matching reference base
        
        int altCount = totalDepth - refDepth; // no. of alternate reads
        
        long double altFreq;
        if (totalDepth == 0) altFreq = 0;
        else altFreq = (long double) altCount / totalDepth;
        
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
        
        // Count stuff idk
        int readSupportCount = 0;
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        if ((rowN+1) % 1000 == 0) printf("Processed row %d\n", rowN+1);
    }
}
