//
//  vcf.hpp
//  MonovarNG
//
//  Created by tianyi on 6/19/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef vcf_hpp
#define vcf_hpp

#include <stdio.h>
#include <vector>
#include <string>

using namespace std;

class VCFDocument {
public:
    void write_stuff(vector<string> bamIDs);
};

#endif /* vcf_hpp */
