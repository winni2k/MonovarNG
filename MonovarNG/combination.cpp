//
//  combination.cpp
//  MonovarNG
//
//  Created by tianyi on 7/4/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "combination.hpp"
#include "wrdouble.hpp"

#include <vector>

using namespace std;

Combination::Combination(){} // default constructor

Combination::Combination(int width): width(width) {
    mem.resize(width+1);
    for (auto &row: mem) row.resize(width+1);
    // first row
    mem[0][0] = 1.0; // 0C0 = 1
    for (int j = 1; j <= width; j++) {
        mem[0][j] = 0.0;
    }
    // other rows
    for (int i = 1; i <= width; i++) {
        mem[i][0] = 1; // nC0 = 1
        for (int j = 1; j <= width; j++) {
            mem[i][j] = mem[i-1][j] + mem[i-1][j-1];
        }
    }
}

const vector<wrdouble> Combination::getRow(int n) const {
    // Gets the row for nC0...nCn
    return mem[n];
}

wrdouble Combination::getValue(int n, int r) const {
    // Gets nCr
    return mem[n][r];
}
