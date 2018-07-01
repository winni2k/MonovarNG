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
    
    
    Config config = setupConfig(argc, argv);
    
    vector<string> bamIDs = getBamIDs(config.bamfileNames);
    
    int numCells = bamIDs.size();
    
    vector<Pileup> pileup = getPileup(numCells, config.pileupFilename);
    
    App app(config, bamIDs, pileup);
    
    app.runAlgo();
    
    return 0;
}
