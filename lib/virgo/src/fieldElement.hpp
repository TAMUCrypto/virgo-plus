#pragma once

#include "constants.h"
#include "typedef.hpp"
#include <cassert>
#include <immintrin.h>
#include <vector>
#include <cstdio>
#include <iostream>

using std::vector;
using std::ostream;

namespace virgo {
#define MASK 4294967295 //2^32-1
#define PRIME 2305843009213693951 //2^61-1

    class fieldElementPacked;

    class fieldElement {
    public:
        fieldElement();

        fieldElement(const fieldElement &b);

        fieldElement(long long x);

        fieldElement(long long x, long long y);

        fieldElement operator+(const fieldElement &other) const;

        fieldElement operator-(const fieldElement &other) const;

        fieldElement operator-() const;

        fieldElement operator*(const fieldElement &other) const;

        bool operator==(const fieldElement &other) const;

        bool operator!=(const fieldElement &other) const;

        fieldElement &operator=(const fieldElement &other);

        fieldElement &operator+=(const fieldElement &other);

        fieldElement &operator-=(const fieldElement &other);

        fieldElement &operator*=(const fieldElement &other);

        friend ostream &operator << (ostream &out, const fieldElement &c);

        bool operator < (const fieldElement &other) const;

        explicit operator bool () const;

        [[nodiscard]] bool isNegative() const;

        [[nodiscard]] unsigned char getBitWidth() const;

        [[nodiscard]] unsigned char getBit(unsigned int i) const;

        [[nodiscard]] __int128_t toint128() const;

        [[nodiscard]] bool isZero();

        [[nodiscard]] fieldElement abs() const;
        [[nodiscard]] fieldElement sqr() const;
        [[nodiscard]] fieldElement inv() const;
        void setAbs();
        void setSqr();
        void setInv();

        void print(::FILE *fileno) const;
        char *toString() const;

        static void init();
        static fieldElement maxWithZero(const fieldElement &a, const fieldElement &b);
        static fieldElement maxUnsigned(const fieldElement &a, const fieldElement &b);
        static fieldElement getRootOfUnity(int log_order); //return a root of unity with log_order 2^[log_order]
        static fieldElement random();

        static fieldElement zero();
        static fieldElement one();
        static vector<fieldElement> generateRandomness(u64 size);
        static fieldElement innerProd(vector<fieldElement>::iterator a, vector<fieldElement>::iterator b, u64 n);

        static fieldElement fastPow(fieldElement x, __uint128_t p);

        static bool initialized;
        static int multCounter, addCounter;
        static bool isCounting;
        static bool isSumchecking;
        static const unsigned long long mod;
        static const int __max_order = 62;

        unsigned long long real;
        unsigned long long img;

    protected:
        static double self_speed_test_mult(int repeat);
        static double self_speed_test_add(int repeat);

        static unsigned long long myMod(unsigned long long x);
        static unsigned long long mymult(const unsigned long long x, const unsigned long long y);
        static unsigned long long randomNumber();

        friend fieldElementPacked;
    };

    class fieldElementPacked {
    public:
        __m256i img, real;

        fieldElementPacked();
        fieldElementPacked(const fieldElement &x0, const fieldElement &x1, const fieldElement &x2, const fieldElement &x3);

        static void init();
        fieldElementPacked operator + (const fieldElementPacked &b) const;
        fieldElementPacked operator * (const fieldElementPacked &b) const;
        fieldElementPacked operator - (const fieldElementPacked &b) const;
        __mmask8 operator == (const fieldElementPacked &b) const;
        __mmask8 operator != (const fieldElementPacked &b) const;
        void getFieldElement(fieldElement *dst) const;

        static const unsigned long long mod;
        static __m256i packed_mod, packed_mod_minus_one;

        static __m256i packed_mymult(const __m256i x, const __m256i y);
        static __m256i packed_myMod(const __m256i x);

    };
}














