//
// Created by 69029 on 6/23/2021.
//

#include "utils.hpp"
#include "config_pc.hpp"

void initHalfTable(vector<F> &beta_f, vector<F> &beta_s, const vector<F>::const_iterator &r, const F &init, u64 first_half, u64 second_half) {
    beta_f.at(0) = init;
    beta_s.at(0) = F_ONE;

    for (u32 i = 0; i < first_half; ++i) {
        for (u32 j = 0; j < (1ULL << i); ++j) {
            auto tmp = beta_f.at(j) * r[i];
            beta_f.at(j | (1ULL << i)) = tmp;
            beta_f.at(j) = beta_f[j] - tmp;
        }
    }

    for (u32 i = 0; i < second_half; ++i) {
        for (u32 j = 0; j < (1ULL << i); ++j) {
            auto tmp = beta_s[j] * r[(i + first_half)];
            beta_s[j | (1ULL << i)] = tmp;
            beta_s[j] = beta_s[j] - tmp;
        }
    }
}

void initBetaTable(vector<F> &beta_g, u8 gLength, const vector<F>::const_iterator &r, const F &init) {
    if (gLength == -1) return;

    static vector<F> beta_f, beta_s;
    int first_half = gLength >> 1, second_half = gLength - first_half;
    u32 mask_fhalf = (1ULL << first_half) - 1;

    assert(beta_g.size() >= 1ULL << gLength);
    myResize(beta_f, 1ULL << first_half);
    myResize(beta_s, 1ULL << second_half);
    if (init != F_ZERO) {
        initHalfTable(beta_f, beta_s, r, init, first_half, second_half);
        for (u32 i = 0; i < (1ULL << gLength); ++i)
            beta_g.at(i) = beta_f.at(i & mask_fhalf) * beta_s.at(i >> first_half);
    } else for (u32 i = 0; i < (1ULL << gLength); ++i)
            beta_g.at(i) = F_ZERO;
}

void initBetaTable(vector<F> &beta_g, int gLength, const vector<F>::const_iterator &r_0,
                   const vector<F>::const_iterator &r_1, const F &alpha, const F &beta) {
    static vector<F> beta_f, beta_s;
    int first_half = gLength >> 1, second_half = gLength - first_half;
    u32 mask_fhalf = (1ULL << first_half) - 1;

    assert(beta_g.size() >= 1ULL << gLength);
    myResize(beta_f, 1ULL << first_half);
    myResize(beta_s, 1ULL << second_half);
    if (beta != F_ZERO) {
        initHalfTable(beta_f, beta_s, r_1, beta, first_half, second_half);
        for (u32 i = 0; i < (1ULL << gLength); ++i)
            beta_g[i] = beta_f[i & mask_fhalf] * beta_s[i >> first_half];
    } else for (u32 i = 0; i < (1ULL << gLength); ++i)
            beta_g[i] = F_ZERO;

    if (alpha == F_ZERO) return;
    initHalfTable(beta_f, beta_s, r_0, alpha, first_half, second_half);
    assert(beta_g.size() >= 1ULL << gLength);
    for (u32 i = 0; i < (1ULL << gLength); ++i)
        beta_g[i] = beta_g[i] + beta_f[i & mask_fhalf] * beta_s[i >> first_half];
}