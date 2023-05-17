//
// Created by 69029 on 7/19/2021.
//

#ifndef VIRGO_PLUS_INPUTCIRCUIT_HPP
#define VIRGO_PLUS_INPUTCIRCUIT_HPP
#include <bits/stdc++.h>
using std::pair;
using std::vector;
using std::map;

typedef long long i64;
typedef unsigned long long u64;
enum gateType {
    Mul, Add, Sub, AntiSub, Naab, AntiNaab, Input, Mulc, Addc, Xor, Not, Copy, SIZE
};

class DAG_gate {
public:
    pair<int, u64> input0, input1;
    bool is_assert;
    gateType ty;
};

#endif //VIRGO_PLUS_INPUTCIRCUIT_HPP
