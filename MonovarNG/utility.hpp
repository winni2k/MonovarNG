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
#include "wrdouble.hpp"

#include <stdio.h>
#include <array>

using namespace std;

namespace utility {
    
    Config setupConfig(int argc, const char *argv[]); // Setup app configurations
    
    vector<string> getBamFilenames(string filename); // Gets bam filenames for all bam files named in file (at filename)
    
    vector<string> getBamIDs(string filename); // Gets bam IDs for all bam files named in file (at filename). RG IDs, not filename
    
    vector<string> readPileup(int numCells, string filename); // Reads from pileup file, saving it as a vector of strings
    
    Pileup getPileup(int numCells, string& row); // Parses a row of pileup and return a pileup object
    
    array<array<array<double, 4>, 4>, 4> genGenotypePriors(double p); // Generates genotype priors matrix given probability p. priors[a][b][c] = p(^ab)(_c)
    
    vector<double> genAltCountPriors(int numCells, double mutationRate = 0.001); // gets the alternate allele count priors, for the given number of cells with reads. Return array size = 2*numCells+1
}
#endif /* utility_hpp */
