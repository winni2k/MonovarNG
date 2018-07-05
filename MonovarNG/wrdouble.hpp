//
//  wrdouble.hpp
//  MonovarNG
//
//  Created by tianyi on 7/3/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#ifndef wrdouble_hpp
#define wrdouble_hpp

#include <stdio.h>
#include <iostream>
#include <cmath>

using namespace std;

struct wrdouble {
    
    // wrdouble = value * base^exp. Can never be negative
    static double base;
    static double invBase; // inverse of base
    
    double value; // 1 <= value < base
    int exponent;
    
    wrdouble(); // default constructor
    wrdouble(double value, int exponent); // constructor given value and exp
    wrdouble(double n); // constructor from double
    
    operator double(); // casting to double
    operator string() const; // casting to string
    
    friend ostream& operator<<(ostream& out, const wrdouble& n); // printing conversion
    
    wrdouble& operator=(double n); // assignment of double
    wrdouble& operator=(const wrdouble& n); // assignment of wrdouble
    
    bool operator<(wrdouble& n); // comparator <
    bool operator>(wrdouble& n); // comparator >
    
    wrdouble operator*(const wrdouble& n) const; // multiplication
    wrdouble operator/(const wrdouble& n) const; // division
    wrdouble operator+(const wrdouble& n) const; // addition
    
    wrdouble operator*(double n) const; // multiplication with double
    wrdouble operator/(double n) const; // division with double
    wrdouble operator+(double n) const; // addition with double
    
    wrdouble& operator*=(const wrdouble& n); // multiplication and assignment
    wrdouble& operator/=(const wrdouble& n); // division and assignment
    wrdouble& operator+=(const wrdouble& n); // addition and assignment
    
    wrdouble& operator*=(double n); // multiplication and assignment with double
    wrdouble& operator/=(double n); // division and assignment with double
    wrdouble& operator+=(double n); // addition and assignment with double
};

#endif /* wrdouble_hpp */
