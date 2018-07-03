//  Wide Range Doubles
//  wrdouble.cpp
//  MonovarNG
//
//  Created by tianyi on 7/3/18.
//  Copyright © 2018 tianyi. All rights reserved.
//

#include "wrdouble.hpp"

#include <cmath>

double wrdouble::base = pow(double(2), 64);
double wrdouble::invBase = pow(double(2), -64);

wrdouble::wrdouble(double value, int exponent): value(value), exponent(exponent) {} // constructor

wrdouble::wrdouble(double n) {
    // constructor from double
    exponent = 0;
    value = n;
    while (value < 1) { // if n < 1
        exponent--;
        value *= base;
    }
    while (value >= base) { // if n >= base
        exponent++;
        value *= invBase; // multiply by invbase to speed up
    }
}


wrdouble::operator double() {
    double n = value;
    n *= pow(base, exponent);
    if (isinf(n)) return 0; // if underflow 
    return n;
}


wrdouble& wrdouble::operator=(double n) {
    // assignment of double
    exponent = 0;
    value = n;
    while (value < 1) { // if n < 1
        exponent--;
        value *= base;
    }
    while (value >= base) { // if n >= base
        exponent++;
        value *= invBase; // multiply by invbase to speed up
    }
    
    return *this;
}

wrdouble& wrdouble::operator=(const wrdouble& n) {
    // assignment of wrdouble
    value = n.value;
    exponent = n.exponent;
    return *this;
}


wrdouble wrdouble::operator*(const wrdouble& n) {
    // multiplication
    int newexp = exponent + n.exponent;
    double newval = value * n.value;
    if (newval >= base) {
        newexp++;
        newval *= invBase;
    }
    return wrdouble(newval, newexp);
}

wrdouble wrdouble::operator/(const wrdouble& n) {
    // division
    int newexp = exponent - n.exponent;
    double newval = value / n.value;
    if (newval < 1) {
        newexp--;
        newval *= base;
    }
    return wrdouble(newval, newexp);
}

wrdouble wrdouble::operator+(const wrdouble& n) {
    // addition
    if (exponent > n.exponent+1) { // if this >> n
        return *this;
    } else if (exponent + 1 < n.exponent) { // if this << n
        return n;
    } else if (exponent == n.exponent) {
        double newval = value + n.value;
        int newexp = exponent;
        if (newval >= base) {
            newexp++;
            newval *= invBase;
        }
        return wrdouble(newval, newexp);
    } else if (exponent == n.exponent + 1) { // this > n, not by too much
        double newval = value + n.value*invBase; // bring n.value to the same exp
        int newexp = exponent;
        if (newval >= base) {
            newexp++;
            newval *= invBase;
        }
        return wrdouble(newval, newexp);
    } else { // this < n, not by too much
        double newval = value*invBase + n.value; // bring value to the same exp
        int newexp = n.exponent;
        if (newval >= base) {
            newexp++;
            newval *= invBase;
        }
        return wrdouble(newval, newexp);
    }
}

wrdouble wrdouble::operator*(double n) {
    // multiplication with double
    return (*this) * wrdouble(n);
}

wrdouble wrdouble::operator/(double n) {
    // division with double
    return (*this) / wrdouble(n);
}

wrdouble wrdouble::operator+(double n) {
    // addition with double
    return (*this) + wrdouble(n);
}


wrdouble& wrdouble::operator*=(const wrdouble& n) {
    // multiplication and assignment
    exponent += n.exponent;
    value *= n.value;
    if (value >= base) {
        exponent++;
        value *= invBase;
    }
    return *this;
}

wrdouble& wrdouble::operator/=(const wrdouble& n) {
    // division and assignment
    exponent -= n.exponent;
    value /= n.value;
    if (value < 1) {
        exponent--;
        value *= base;
    }
    return *this;
}

wrdouble& wrdouble::operator+=(const wrdouble& n) {
    // addition and assignment
    if (exponent > n.exponent+1) {} // if this >> n
    else if (exponent + 1 < n.exponent) { // if this << n
        value = n.value;
        exponent = n.exponent;
    } else if (exponent == n.exponent) {
        value += n.value;
        if (value >= base) {
            exponent++;
            value *= invBase;
        }
    } else if (exponent == n.exponent + 1) { // this > n, not by too much
        value += n.value*invBase; // bring n.value to the same exp
        if (value >= base) {
            exponent++;
            value *= invBase;
        }
    } else { // this < n, not by too much
        value = value*invBase + n.value; // bring value to the same exp
        exponent = n.exponent;
        if (value >= base) {
            exponent++;
            value *= invBase;
        }
    }
    
    return *this;
}

wrdouble& wrdouble::operator*=(double n) {
    // multiplication and assignment with double
    return (*this) *= wrdouble(n);
}

wrdouble& wrdouble::operator/=(double n) {
    // division and assignment with double
    return (*this) /= wrdouble(n);
}

wrdouble& wrdouble::operator+=(double n) {
    // addition and assignment with double
    return (*this) += wrdouble(n);
}
