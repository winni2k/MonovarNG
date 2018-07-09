//
//  phred.cpp
//  MonovarNG
//
//  Created by tianyi on 7/6/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "phred.hpp"

#include <cmath>

using namespace std;

Phred::Phred() {
    // Default constructor
    for (int i = 0; i <= 250; i++) {
        qualities[i] = pow(double(10), double(-i)/10);
    }
}
