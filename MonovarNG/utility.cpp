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

Config utility::setupConfig(int argc, const char *argv[]) {
    // Setup app configurations
    Config config;
    
    // basic implementation, only parsing key variables
    if (argc != 5) {
        throw invalid_argument("Incorrect number of arguments; received " + to_string(argc-1) + " instead of 4.");
    }
    
    config.referenceFilename = argv[1];
    config.bamfileNames = argv[2];
    config.pileupFilename = argv[3];
    config.outputFilename = argv[4];
    
    return config;
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

vector<Pileup> utility::getPileup(int numCells, string filename) {
    // Gets a list of pileup rows
    ifstream pileupFile;
    pileupFile.open(filename);
    
    vector<Pileup> rows;
    
    string row;
    while (getline(pileupFile, row)) {
        boost::trim(row);
        if (row.size()) {
            rows.push_back(Pileup(numCells, row));
        }
    }
    
//    cout << to_string(rows.size()) << " rows read." << endl;
    
    return rows;
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

