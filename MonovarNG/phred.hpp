//
//  phred.hpp
//  MonovarNG
//
//  Created by tianyi on 7/6/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef phred_hpp
#define phred_hpp

#include <stdio.h>

struct Phred {
    double qualities[300];
    Phred(); // default constructor
};

#endif /* phred_hpp */
