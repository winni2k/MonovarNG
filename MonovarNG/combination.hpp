//
//  combination.hpp
//  MonovarNG
//
//  Created by tianyi on 7/4/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef combination_hpp
#define combination_hpp

#include "wrdouble.hpp"

#include <stdio.h>
#include <vector>

using namespace std;

class Combination { // Computes nCr and stores it for use
    int width;
    vector<vector<wrdouble>> mem; // stores nCr
public:
    Combination(); // default constructor
    Combination (int width); // width = largest n, 0C0 to widthCwidth. Width at least 1
    const vector<wrdouble> getRow(int n) const; // gets the row for nC0...nCn
    wrdouble getValue(int n, int r) const; // gets nCr
};

#endif /* combination_hpp */
