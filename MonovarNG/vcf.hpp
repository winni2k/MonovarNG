//
//  vcf.hpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef vcf_hpp
#define vcf_hpp

#include "wrdouble.hpp"

#include <stdio.h>
#include <vector>
#include <string>
#include <array>
#include <fstream>

using namespace std;

class VCFDocument {
private:
    ofstream outputFile;
public:
    VCFDocument(string filename); // Initialization function, sets up output file
    void writeDefHeader(); // writes default header of vcf file, containing date and format specs
    void writeHeaderInfo(string referenceFilename, vector<string> bamIDs); // writes specific info, like reference file, column headers
    void writeRow(string chromosome, int posID, char ref, char alt, double quality, double wilcoxon, double qualityByDepth, double strandBias, double psarr, vector<int> genotypes, int depth, vector<pair<int, int>> cellDepths, vector<array<wrdouble, 3>> likelihoods); // writes a row, for mutation at a given site into file
};

#endif /* vcf_hpp */
