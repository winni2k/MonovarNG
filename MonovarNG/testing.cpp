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
#include "wrdouble.hpp"

#include <iostream>
#include <array>
#include <chrono>
#include <vector>

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
    
//    wrdouble a = 1e-100;
//    a *= 3e50;
//    cout << a;
//    
//    wrdouble a = 0;
//    wrdouble b = wrdouble(1.1, -50);
//    cout << b << endl;
//    cout << a+b << endl;
//    cout << wrdouble(0) + b << endl;
    // testing double multiplication
    
//    vector<double> a;
//    for (int i = 0; i < 1000000; i++) a.push_back(1.0);
//    
//    auto start = chrono::high_resolution_clock::now();
//    for (int i = 0; i < 1000000; i++) {
//        for (int j = 0; j < 1000; j++) {
//            a[i] *= 1.23456789;
//        }
//    }
//    auto end = chrono::high_resolution_clock::now();
//    auto elapsed = end-start;
//    printf("Time for 1e9 double multiplications = %lld\n", chrono::duration_cast<chrono::milliseconds>(elapsed).count());
//    
//    // testing wrdouble multiplication
//    vector<wrdouble> b;
//    for (int i = 0; i < 1000000; i++) b.push_back(wrdouble(1.0));
//    wrdouble mult = 1.23456789;
//    
//    start = chrono::high_resolution_clock::now();
//    for (int i = 0; i < 1000000; i++) {
//        for (int j = 0; j < 1000; j++) {
//            b[i] = b[i] * mult;
//        }
//    }
//    end = chrono::high_resolution_clock::now();
//    elapsed = end-start;
//    printf("Time for 1e9 wrdouble multiplications = %lld\n", chrono::duration_cast<chrono::milliseconds>(elapsed).count());
//    
//    // testing double addition
//    
//    for (int i = 0; i < 1000000; i++) a[i] = 1.0;
//    
//    start = chrono::high_resolution_clock::now();
//    for (int i = 0; i < 1000000; i++) {
//        for (int j = 0; j < 1000; j++) {
//            a[i] += 1.23456789;
//        }
//    }
//    end = chrono::high_resolution_clock::now();
//    elapsed = end-start;
//    printf("Time for 1e9 double additions = %lld\n", chrono::duration_cast<chrono::milliseconds>(elapsed).count());
//    
//    // testing wrdouble addition
//    for (int i = 0; i < 1000000; i++) b[i] = 1.0;
//    wrdouble add = 1.23456789;
//    
//    start = chrono::high_resolution_clock::now();
//    for (int i = 0; i < 1000000; i++) {
//        for (int j = 0; j < 1000; j++) {
//            b[i] = b[i]+add;
//        }
//    }
//    end = chrono::high_resolution_clock::now();
//    elapsed = end-start;
//    printf("Time for 1e9 wrdouble additions = %lld\n", chrono::duration_cast<chrono::milliseconds>(elapsed).count());
    
    // testing wrdouble assignment
    vector<wrdouble> b;
    
    for (int i = 0; i < 1000000; i++) b.push_back(1.0);
    wrdouble add = 1.23456789;
    
    auto start = chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000000; i++) {
        for (int j = 0; j < 1000; j++) {
            b[i] = add;
        }
    }
    auto end = chrono::high_resolution_clock::now();
    auto elapsed = end-start;
    printf("Time for 1e9 wrdouble assignments = %lld\n", chrono::duration_cast<chrono::milliseconds>(elapsed).count());

}
