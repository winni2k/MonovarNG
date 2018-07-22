//
//  utility.cpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "utility.hpp"
#include "pileup.hpp"
#include "wrdouble.hpp"

#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>

using namespace std;
using namespace utility;

Config utility::setupConfig(int argc, const char *argv[]) {
    // Setup app configurations
    Config config;
    
    if (argc < 5) {
        throw invalid_argument("Incorrect arguments.\nUsage: monovar referenceFile bamFilenames pileupFile outputFile [-patm]\nOptions:\n-t: Threshold to be used for variant calling (Recommended value: 0.05)\n-p: Offset for prior probability for false-positive error (Recommended value: 0.002)\n-a: Offset for prior probability for allelic drop out (Default value: 0.2)\n-m: Number of threads to use in multiprocessing (Default value: 4)");
    }
    
    config.referenceFilename = argv[1];
    config.bamfileNames = argv[2];
    config.pileupFilename = argv[3];
    config.outputFilename = argv[4];
    
    for (int i = 5; i < argc; i++) {
        string token = argv[i];
        if (token.size() > 0 && token[0] == '-') {
            if (token.size() < 2) continue; // no option given
            if (i+1 == argc) continue; // no argument given
            switch(token[1]) {
                case 't':
                    config.mutationThreshold = atof(argv[i+1]);
                    break;
                case 'p':
                    config.pFalsePositive = atof(argv[i+1]);
                    break;
                case 'a':
                    config.pDropout = atof(argv[i+1]);
                    break;
                case 'm':
                    config.numThreads = atoi(argv[i+1]);
                    break;
            }
        }
    }
    
    return config;
}

vector<string> utility::getBamFilenames(string filename) {
    // Gets bam filenames for all bam files named in file (at filename)
    ifstream bamNameFile;
    bamNameFile.open(filename);
    
    // Get bam names
    vector<string> bamNames;
    string tmp;
    while (getline(bamNameFile, tmp)) {
        boost::trim(tmp);
        if (tmp.size()) { // check if tmp is whitespace only 
            bamNames.push_back(tmp);
        }
    }
    bamNameFile.close();
    return bamNames;
}

vector<string> utility::getBamIDs(string filename) {
    // Gets bam IDs for all bam files named in file (at filename)
    vector<string> ids;
    
    ifstream bamNameFile;
    bamNameFile.open(filename);
    
    // Get bam names
    vector<string> bamNames;
    string tmp;
    while (getline(bamNameFile, tmp)) {
        boost::trim(tmp);
        if (tmp.size()) { // check if tmp is whitespace only 
            bamNames.push_back(tmp);
        }
    }
    bamNameFile.close();

    // Open each bam file and get its readgroup ID
    for (string file: bamNames) {
        samFile * samfile = sam_open(file.c_str(), "r");
        bam_hdr_t * header = sam_hdr_read(samfile);
        string headerText = (*header).text;
        
        vector<string> lines;
        boost::split(lines, headerText, boost::is_any_of("\n"));
        for (string line: lines) {
            if (boost::starts_with(line, "@RG")) {
                vector<string> tokens;
                boost::split(tokens, line, boost::is_any_of("\t"));
                for (string token: tokens) {
                    if (boost::starts_with(token, "ID")) {
                        vector<string> items;
                        boost::split(items, token, boost::is_any_of(":"));
                        string id = items[1];
                        boost::trim(id);
                        ids.push_back(id);
                    }
                }
            }
        }
    }
    
    return ids;
}

vector<string> utility::readPileup(int numCells, string filename) {
    // Reads from pileup file, saving it as a vector of strings
    printf("Reading from %s\n", filename.c_str());
    ifstream pileupFile;
    pileupFile.open(filename);
    
    vector<string> rows;
    
    string row;
    while (getline(pileupFile, row)) {
        boost::trim(row);
        if (row.size()) {
            rows.push_back(row);
        }
    }
    printf("%d positions read.\n", rows.size());
    return rows;
}

Pileup utility::getPileup(int numCells, string& row) {
    // Parses a row of pileup and return a pileup object
    return Pileup(numCells, row);
}

array<array<array<double, 4>, 4>, 4> utility::genGenotypePriors(double p) {
    // Generates genotype priors matrix given probability p. priors[a][b][c] = p(^ab)(_c)
    array<array<array<double, 4>, 4>, 4> priors;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                if (i == j) {
                    if (i == k) { // AA A, CC C...
                        priors[i][j][k] = 1 - 3*p;
                    } else { // AA C, AA T...
                        priors[i][j][k] = p;
                    }
                } else {
                    if (i == k || j == k) { // AC A, AC C...
                        priors[i][j][k] = (1-2*p)/2;
                    } else { // AC T, ACG...
                        priors[i][j][k] = p;
                    }
                }
            }
        }
    }
    return priors;
}

vector<double> utility::genAltCountPriors(int numCells, double mutationRate) {
    // gets the alternate allele count priors, for the given number of cells with reads. Return array size = 2*numCells+1
    vector<double> priors; 
    for (int l = 0; l <= 2*numCells; l++) {
        double prob = 0;
        for (int i = 1; i <= 2*numCells-1; i++) {
            prob += double(1)/i;
        }
        prob = (1-mutationRate*prob)/2;
        
        if (l == 0 || l == 2*numCells) priors.push_back(prob);
        else priors.push_back(mutationRate/l);
    }
    
    return priors;
}

