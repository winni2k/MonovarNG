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

wrdouble::wrdouble() {} // default constructor

wrdouble::wrdouble(double value, int exponent): value(value), exponent(exponent) {} // constructor

wrdouble::wrdouble(double n) {
    // constructor from double
    if (n == 0) { // if zero
        value = 0;
        exponent = -1e6; // this number is very small!
        return;
    }
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
    // casting to double
    double n = value;
    n *= pow(base, exponent);
    if (isinf(n) || isnan(n)) return 0; // if underflow 
    return n;
}

wrdouble::operator string() const {
    // casting to string
    if (value == 0) return "0";
    double logged = log10(value) + log10(2) * 64 * exponent;
    double val = pow(double(10), logged-floor(logged));
    int exp = floor(logged);
    if (abs(exp) <= 5) return to_string(val*pow(10, exp));
    else return to_string(val) + "e" + to_string(exp);
}

ostream& operator<<(ostream& out, const wrdouble& n) {
    // printing conversion
    out << string(n); 
    return out;
}

wrdouble& wrdouble::operator=(double n) {
    // assignment of double
    if (n == 0) { // if zero
        value = 0;
        exponent = -1e6; // this number is very small!
        return *this;
    }
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


bool wrdouble::operator<(wrdouble& n) {
    // comparator <
    if (exponent < n.exponent) return true;
    else if (exponent > n.exponent) return false;
    else return value < n.value;
}

bool wrdouble::operator>(wrdouble& n) {
    // comparator >
    if (exponent > n.exponent) return true;
    else if (exponent < n.exponent) return false;
    else return value > n.value;
}


wrdouble wrdouble::operator*(const wrdouble& n) const {
    // multiplication
    int newexp = exponent + n.exponent;
    double newval = value * n.value;
    if (newval >= base) {
        newexp++;
        newval *= invBase;
    }
    return wrdouble(newval, newexp);
}

wrdouble wrdouble::operator/(const wrdouble& n) const {
    // division
    int newexp = exponent - n.exponent;
    double newval = value / n.value;
    if (newval < 1 && newval != 0) {
        newexp--;
        newval *= base;
    }
    return wrdouble(newval, newexp);
}

wrdouble wrdouble::operator+(const wrdouble& n) const {
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

wrdouble wrdouble::operator*(double n) const {
    // multiplication with double
    return (*this) * wrdouble(n);
}

wrdouble wrdouble::operator/(double n) const {
    // division with double
    return (*this) / wrdouble(n);
}

wrdouble wrdouble::operator+(double n) const {
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
    if (value < 1 && value != 0) {
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
