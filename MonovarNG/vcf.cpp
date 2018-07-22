//
//  vcf.cpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "vcf.hpp"
#include "wrdouble.hpp"

#include <stdio.h>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>

using namespace std;

VCFDocument::VCFDocument(string filename) {
    // Initialization function, sets up output file
    outputFile.open(filename);
}

void VCFDocument::writeDefHeader() {
    // writes header of vcf file
    outputFile << "##fileformat=VCFv4.1" << endl;
    auto time = chrono::system_clock::now();
    time_t time_c = chrono::system_clock::to_time_t(time);
    outputFile << "##fileDate=" << put_time(localtime(&time_c), "%F") << endl;
    outputFile << "##source=MonoVar" << endl;
    outputFile << "##FILTER=<ID=LowQual,Description=\"Low quality\">" << endl;
    outputFile << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">" << endl;
    outputFile << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">" << endl;
    outputFile << "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">" << endl;
    outputFile << "##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">" << endl;
    outputFile << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">" << endl;
    outputFile << "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">" << endl;
    outputFile << "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">" << endl;
    outputFile << "##INFO=<ID=MPR,Number=1,Type=Float,Description=\"Log Odds Ratio of maximum value of probability of observing non-ref allele to the probability of observing zero non-ref allele\">" << endl;
    outputFile << "##INFO=<ID=PSARR,Number=1,Type=Float,Description=\"Ratio of per-sample Alt allele supporting reads to Ref allele supporting reads\">" << endl;
    outputFile << "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << endl;
    outputFile << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">" << endl;
    outputFile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">" << endl;
    outputFile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
    outputFile << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">" << endl;
}

void VCFDocument::writeHeaderInfo(string referenceFilename, vector<string> bamIDs) {
    outputFile << "##reference=file:" << referenceFilename << endl;
    outputFile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
    for (auto bam: bamIDs) {
        outputFile << "\t" << bam;
    }
    outputFile << endl;
}

void VCFDocument::writeRow(string chromosome, int posID, char ref, char alt, double quality, double wilcoxon, double qualityByDepth, double strandBias, double psarr, vector<int> genotypes, int depth, vector<pair<int, int>> cellDepths, vector<array<wrdouble, 3>> likelihoods) {
    // writes a row, for mutation at a given site into file
    char baseMap[5] = {'A', 'C', 'T', 'G'};
    outputFile << chromosome << "\t" << posID << "\t.\t" << baseMap[ref] << "\t" << baseMap[alt] << "\t" << quality << "\t.\t"; 
    
    // compute alt count and stuff
    int altCount = 0, alleleCount = 0;
    for (int i: genotypes) {
        if (i != -1) {
            altCount += i;
            alleleCount += 2;
        }
    }
    double altFreq = double(altCount)/alleleCount;
    outputFile << "AC=" << altCount << ";AF=" << altFreq << ";AN=" << alleleCount << ";";
    
    // wilcoxon
    outputFile << "BaseQRankSum=" << wilcoxon << ";";
    
    // depth
    outputFile << "DP=" << depth << ";";
    
    // quality by depth
    outputFile << "QD=" << qualityByDepth << ";";
    
    // strand bias
    outputFile << "SOR=" << strandBias << ";";
    
    
    // psarr
    outputFile << "PSARR=" << psarr;
    
    // Individual cell format
    outputFile << "\tGT:AD:DP:GQ:PL";
    int likelihoodsIndex = 0;
    for (int i = 0; i < genotypes.size(); i++) {
        if (genotypes[i] == -1) outputFile << "\t./.";
        else {
            outputFile << "\t";
            if (genotypes[i] == 0) outputFile << "0/0";
            else if (genotypes[i] == 1) outputFile << "0/1";
            else if (genotypes[i] == 2) outputFile << "1/1";
            
            outputFile << ":" << cellDepths[i].first << "," << cellDepths[i].second;
            outputFile << ":" << cellDepths[i].first+cellDepths[i].second;
            
            double quals[3];
            for (int j = 0; j < 3; j++) quals[j] = likelihoods[likelihoodsIndex][j].phred();
            double lowest = 1e9;
            for (double qual: quals) lowest = min(lowest, qual);
            for (int j = 0; j < 3; j++) quals[j] = round(quals[j]-lowest);
            int secondLowest = 1e9;
            for (int j = 0; j < 3; j++) if (quals[j]) secondLowest = min(secondLowest, int(quals[j]));
            
            outputFile << ":" << secondLowest << ":";
            for (int j = 0; j < 3; j++) {
                if (j != 0) outputFile << ",";
                outputFile << quals[j];
            }
            
            likelihoodsIndex++;
        }
    }
    
    // Final genotype summary
    outputFile << "\t<";
    for (int i: genotypes) {
        if (i == -1) outputFile << "X";
        else outputFile << i;
    }
    outputFile << ">";
    outputFile << endl;
}
