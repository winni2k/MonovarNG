//
//  main.cpp
//  MonovarNG
//
//  Created by tianyi on 6/9/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "testing.hpp"
#include "utility.hpp"
#include "config.hpp"
#include "app.hpp"
#include "pileup.hpp"

#include <string>
#include <vector>

using namespace std;
using namespace utility;

int main(int argc, const char * argv[]) {
//    test();
    
    auto start = chrono::high_resolution_clock::now();
    Config config = setupConfig(argc, argv);
    
    vector<string> bamIDs = getBamIDs(config.bamfileNames);
    
    int numCells = bamIDs.size();
    
    vector<string> pileup = readPileup(numCells, config.pileupFilename);
    
    App app(config, bamIDs, pileup);
    
    auto end = chrono::high_resolution_clock::now();
    auto setupTime = end-start;
    printf("Time for setup = %lldms\n", chrono::duration_cast<chrono::milliseconds>(setupTime).count());
    
    start = chrono::high_resolution_clock::now();
    app.runAlgo();
    end = chrono::high_resolution_clock::now();
    auto algoTime = end-start;
    printf("Time for algo = %lldms\n", chrono::duration_cast<chrono::milliseconds>(algoTime).count());
    printf("Total time = %lfs\n", double(chrono::duration_cast<chrono::milliseconds>(setupTime+algoTime).count())/1000);
    
    exit(0); // skip deallocation
    return 0;
}
