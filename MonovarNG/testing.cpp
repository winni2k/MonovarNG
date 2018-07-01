//
//  testing.cpp
//  MonovarNG
//
//  Created by tianyi on 6/27/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "testing.hpp"
#include "single_cell_pos.hpp"
#include "utility.hpp"

#include <iostream>
#include <array>

using namespace std;
using namespace utility;

void test() {
    // Universal testing function
    
    // Cell testing
//    auto cell = SingleCellPos(205, "......C.................................................................................................C....................................................................................................", "HFFHHHHHH2HHFHFGFHFHHFHHHHHHHHH5HGGF@EHHHHHF2HHHFHHH5BHHHHHHHGHHHHHHHHHHHDHHHGGFHEHHHHHHHH@HHHHHHHHHHHHHHGGGHFH@HHHHHHHH5HHHFGHDHHFHHGGDHG5FF5FHHHGFBFHHFFHHHHHHHHG5HHHBHGHBHHDHHHHHHH2HHFHHHHHHHHHGHHDFFHGGH");
//    cell.removeInsDels();
//    cell.removeStartEnd();
//    cell.cleanupBases('G');
//    cell.truncateReads();
//    cell.computeQuality();
//    
//    array<int, 4> baseFreq = cell.baseFreq();
//    for (int i = 0; i < 4; i++) cout << baseFreq[i] << endl;
////    cout << cell.numReads;
////    for (int i = 0; i < cell.qualityString.size(); i++) {
////        printf("%c - %c - %d - %lf\n", cell.bases[i], cell.qualityString[i], (int)cell.qualityString[i], cell.qualities[i]);
////    }
    
//    auto matrix = genGenotypePriors(0.1);
//    string m = "ACTG";
//    for (int i = 0; i < 4; i++) {
//        for (int j = 0; j < 4; j++) {
//            printf("%c%c -", m[i], m[j]);
//            for (int k = 0; k < 4; k++) {
//                printf(" %c = %lf", m[k], matrix[i][j][k]);
//            }
//            printf("\n");
//        }
//    }
}
