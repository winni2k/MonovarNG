//
//  config.hpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef config_hpp
#define config_hpp

#include <stdio.h>
#include <string>

struct Config {
    // Main configuration for CLI
    std::string referenceFilename; // reference filename
    std::string bamfileNames; // name of file containing bamfile names
    std::string pileupFilename; // name of pileup file
    std::string outputFilename; // name of output file
    
    double mutationThreshold = 0.05; // threshold for variant calling
    double pFalsePositive = 0.002; // p_e, prior probability for false positive 
    double pDropout = 0.02; // p_ad, prior probability for allelic dropout
    
    int numThreads = 1; // number of threads for multiprocessing
    
    bool useConsensusFilter = false; // whether to use Consensus Filter (CF) 
};

#endif /* config_hpp */
