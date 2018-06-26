//
//  utility.cpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "utility.hpp"
#include "pileup.hpp"

#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
