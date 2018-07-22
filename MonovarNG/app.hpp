//
//  app.hpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef app_hpp
#define app_hpp

#include "config.hpp"
#include "vcf.hpp"
#include "pileup.hpp"
#include "utility.hpp"
#include "combination.hpp"

#include <stdio.h>
#include <mutex>

using namespace std;
using namespace utility;

class App {
    // Main application, controls the algorithm flow
    
    double mutationThreshold; // threshold for variant calling
    double pFalsePositive; // p_e, prior probability for false positive 
    double pDropout; // p_ad, prior probability for allelic dropout
    
    int numThreads; // number of threads for multiprocessing
    
    bool useConsensusFilter; // whether to use Consensus Filter (CF) 
    
    int numCells; // number of cells processed
    int numPos; // number of positions being processed
    
    VCFDocument output;
    mutex outputMutex;
    
    Combination combi; // computes nCr
    Phred phred; // computes phred probabilities
    
    vector<string>& pileup;
//    vector<Pileup>& positions;
    
public:
    App(Config& config, vector<string>& bamIDs, vector<string>& pileupRows);
    
    void processRow(int rowN); // processes row of data
    void runAlgo(); // Runs main algorithm
};

#endif /* app_hpp */
