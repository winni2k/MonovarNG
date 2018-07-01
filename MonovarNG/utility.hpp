//
//  utility.hpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef utility_hpp
#define utility_hpp

#include "config.hpp"
#include "pileup.hpp"

#include <stdio.h>
#include <array>

using namespace std;

namespace utility {
    
    Config setupConfig(int argc, const char *argv[]); // Setup app configurations
    
    vector<string> getBamIDs(string filename); // Gets bam IDs for all bam files named in file (at filename)
    
    vector<Pileup> getPileup(int numCells, string filename); // Gets a list of pileup rows
    
    array<array<array<double, 4>, 4>, 4> genPriorMatrix(double p); // Generates priors matrix given probability p. priors[a][b][c] = p(^ab)(_c)
}
#endif /* utility_hpp */
