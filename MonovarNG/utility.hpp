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
    
    array<array<array<double, 4>, 4>, 4> genGenotypePriors(double p); // Generates genotype priors matrix given probability p. priors[a][b][c] = p(^ab)(_c)
    
    vector<double> genAltCountPriors(int numCells, double mutationRate = 0.001); // gets the alternate allele count priors, for the given number of cells with reads. Return array size = 2*numCells+1
}
#endif /* utility_hpp */
