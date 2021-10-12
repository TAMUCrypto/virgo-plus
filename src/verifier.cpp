#include "verifier.h"
#include "utils.hpp"
#include <bits/stdc++.h>
using std::unique_ptr;
using std::make_unique;
using std::fill;
using std::cerr;
using std::endl;
using std::max;
using std::min;

verifier::verifier(prover *pr, const layeredCircuit &cir) :
        p(pr), C(cir) {
    final_claims_v.resize(C.size);
    for (int i = 1; i < C.size; ++i)
        final_claims_v[i].resize(i);

    for (auto & i : coeff_r) i.resize(C.size);

    r_v.resize(C.size + 2);

    // make the prover ready
    assert(p != NULL);
    p->init();

    int max_bl = 0;
    for (auto &c: C.circuit) max_bl = max(max_bl, c.bitLength);
    int max_dad_bl = 0;
    for (auto &c: C.circuit) max_dad_bl = max(max_dad_bl, c.maxDadBitLength);
    beta_g.resize(1ULL << max(max_bl, max_dad_bl));
    beta_u.resize(1ULL << max_bl);
    beta_v.resize(1ULL << max_bl);

    r_u.resize(max_bl);
    r_liu.resize(max_bl);
    for (int i = 1; i < C.size; ++i) {
        if (~C.circuit[i].maxDadBitLength)
            r_v[i].resize(C.circuit[i].maxDadBitLength);
    }
    sig.resize(C.size);

    initBetaTable(beta_g, max(max_bl, max_dad_bl), r_u.begin(), F::zero());

#ifdef USE_VIRGO
    commit_pt = commit_vt = 0.0;
    commit_ps = 0;
#endif
}

void verifier::betaInitPhase1(int depth, const F &assert_random) {
    int cur_bl = C.circuit[depth].bitLength, pre_bl = C.circuit[depth - 1].bitLength;
    initBetaTable(beta_g, cur_bl, r_liu.begin(), F_ONE);
    for (i64 g = 0; g < C.circuit[depth].size; ++g)
        if (C.circuit[depth].gates[g].is_assert) beta_g[g] *= assert_random;
    initBetaTable(beta_u, pre_bl, r_u.begin(), F_ONE);
}

void verifier::betaInitPhase2(int depth) {
    int pre_bl = C.circuit[depth].maxDadBitLength;
    initBetaTable(beta_v, pre_bl, r_v[depth].begin(), F_ONE);
}

void verifier::predicatePhase1(int layer_id) {
    auto &cur = C.circuit[layer_id];
    coeff_l[(u64) gateType::Copy] = F_ZERO;
    coeff_l[(u64) gateType::Not] = F_ZERO;
    coeff_l[(u64) gateType::Addc] = F_ZERO;
    coeff_l[(u64) gateType::Mulc] = F_ZERO;
    bias = F_ZERO;
    for (u64 g = 0; g < cur.size; ++g) {
        auto &gate = cur.gates[g];
        switch (gate.ty) {
            case gateType::Addc:
                bias += beta_g[g] * beta_u[gate.u] * gate.c;
            case gateType::Not: case gateType::Copy:
                coeff_l[(u64) gate.ty] += beta_g[g] * beta_u[gate.u];
                break;
            case gateType::Mulc:
                coeff_l[(u64) gate.ty] += beta_g[g] * beta_u[gate.u] * gate.c;
        }
    }

    fill(coeff_r[(u64) gateType::Add].begin(), coeff_r[(u64) gateType::Add].end(), F_ZERO);
    fill(coeff_r[(u64) gateType::Sub].begin(), coeff_r[(u64) gateType::Sub].end(), F_ZERO);
    fill(coeff_r[(u64) gateType::AntiSub].begin(), coeff_r[(u64) gateType::AntiSub].end(), F_ZERO);
    fill(coeff_r[(u64) gateType::Mul].begin(), coeff_r[(u64) gateType::Mul].end(), F_ZERO);
    fill(coeff_r[(u64) gateType::Naab].begin(), coeff_r[(u64) gateType::Naab].end(), F_ZERO);
    fill(coeff_r[(u64) gateType::AntiNaab].begin(), coeff_r[(u64) gateType::AntiNaab].end(), F_ZERO);
    fill(coeff_r[(u64) gateType::Xor].begin(), coeff_r[(u64) gateType::Xor].end(), F_ZERO);
}

void verifier::predicatePhase2(int layer_id) {
    auto &cur = C.circuit[layer_id];

    coeff_l[(u64) gateType::Copy] *= beta_v[0];
    coeff_l[(u64) gateType::Not] *= beta_v[0];
    coeff_l[(u64) gateType::Addc] *= beta_v[0];
    coeff_l[(u64) gateType::Mulc] *= beta_v[0];
    bias *= beta_v[0];
    for (u64 g = 0; g < cur.size; ++g) {
        auto &gate = cur.gates[g];
        switch (gate.ty) {
            case gateType::Add:
            case gateType::Sub:
            case gateType::AntiSub:
            case gateType::Mul:
            case gateType::Naab:
            case gateType::AntiNaab:
            case gateType::Xor:
                coeff_r[(u64) gate.ty][gate.l] += beta_g[g] * beta_u[gate.u] * beta_v[gate.lv];
        }
    }
}

F verifier::getFinalValue(int layer_id, const F &claim_u,
                                     vector<F>::const_iterator claim_v) {
    auto res = coeff_l[(u64) gateType::Not] * (F_ONE - claim_u)
               + coeff_l[(u64) gateType::Copy] * claim_u
               + coeff_l[(u64) gateType::Addc] * claim_u + bias
               + coeff_l[(u64) gateType::Mulc] * claim_u;
    for (int j = 0; j < layer_id; ++j) {
        auto tmp = coeff_r[(u64) gateType::Add][j] * (claim_u + claim_v[j])
                   + coeff_r[(u64) gateType::Sub][j] * (claim_u - claim_v[j])
                   + coeff_r[(u64) gateType::AntiSub][j] * (claim_v[j] - claim_u)
                   + coeff_r[(u64) gateType::Mul][j] * (claim_u * claim_v[j])
                   + coeff_r[(u64) gateType::Naab][j] * (claim_v[j] - claim_u * claim_v[j])
                   + coeff_r[(u64) gateType::AntiNaab][j] * (claim_u - claim_u * claim_v[j])
                   + coeff_r[(u64) gateType::Xor][j] * (claim_u + claim_v[j] - F(2) * claim_u * claim_v[j]);
        res = res + tmp;
    }
    return res;
}

bool verifier::verify()
{
#ifdef USE_VIRGO
    auto merkle_root_l = p->commit_private();
    poly_ver.p = &(p->poly_prover);
#endif

    verify_timer.start();
    verify_slow_timer.start();

    for (int i = 0; i < C.circuit[C.size - 1].bitLength; ++i)
        r_liu[i] = F::random();
    vector<F>::const_iterator r_0 = r_liu.begin();

    verify_timer.stop();
    verify_slow_timer.stop();

	F previousSum = p->Vres(r_0, C.circuit[C.size - 1].bitLength);
	p -> sumcheckInitAll(r_0);

	for (int i = C.size - 1; i; --i)
	{
        p->sumcheckInit();
        if (!verifyPhase1(i, previousSum)) return false;
        if ((~C.circuit[i].maxDadBitLength) && !verifyPhase2(i, previousSum)) return false;

        auto test_value = getFinalValue(i, final_claim_u, final_claims_v[i].begin());
		if (previousSum != test_value)
		{
            fprintf(stderr, "test: %s previous:%s\n", test_value.toString(), previousSum.toString());
			fprintf(stderr, "Verification fail, semi final, circuit level %d\n", i);
			return false;
		}

        if (!verifyLiu(i, previousSum)) return false;
	}

//    Check the correctness of the input by open polynomial commitment on a random point.
#ifdef USE_VIRGO
	if (!verifyPoly(merkle_root_l, previousSum)) return false;
#endif

    fprintf(stderr, "Verification pass\n");
    fprintf(stdout, "Input size %d\n", C.circuit[0].size);
    fprintf(stdout, "Prove Time %lf\n", p->proveTime());

    fprintf(stdout, "verify time %lf = %lf + %lf(slow)\n", verifySlowTime() + verifyTime(), verifyTime(), verifySlowTime());
    fprintf(stdout, "proof size = %lf kb\n", p->proofSize());
#ifdef USE_VIRGO
    fprintf(stdout, "Polynomial commitment: prove time %lf, verify time %lf, proof size %lf kb\n", p->poly_prover.total_time + commit_pt, commit_vt, commit_ps / 1024.0);
#endif
#ifdef USE_HYRAX_P224
    fprintf(stdout, "Polynomial commitment: prove time %lf, verify time %lf, proof size %lf kb\n", p -> polyProverTime(), poly_v -> getVT(), p -> polyProofSize());
#endif
    return true;
}

bool verifier::verifyPhase1(int layer_id, F &previousSum) {
    verify_timer.start();
    verify_slow_timer.start();

    auto &cur = C.circuit[layer_id], &pre = C.circuit[layer_id - 1];
    for (auto &x: r_u) x = F::random();
    F previousRandom = F_ZERO;

    verify_timer.stop();
    verify_slow_timer.stop();

    F assert_random = F::random();
    p->sumcheckInitPhase1(assert_random);
    for (int j = 0; j < pre.bitLength; ++j)
    {
        auto poly = p->sumcheckUpdatePhase1(previousRandom);
        verify_timer.start();
        verify_slow_timer.start();
        if (poly.eval(0) + poly.eval(1) != previousSum)
        {
            fprintf(stderr, "Verification fail, phase1, circuit %d, current bit %d\n", layer_id, j);
            return false;
        }
        previousRandom = r_u[j];
        previousSum = poly.eval(r_u[j]);
        verify_timer.stop();
        verify_slow_timer.stop();
    }

    p ->sumcheckFinalize1(previousRandom, final_claim_u);

    verify_slow_timer.start();

    betaInitPhase1(layer_id, assert_random);
    predicatePhase1(layer_id);

    verify_slow_timer.stop();
    return true;
}

bool verifier::verifyPhase2(int layer_id, F &previousSum) {
    verify_timer.start();
    verify_slow_timer.start();

    auto &cur = C.circuit[layer_id], &pre = C.circuit[layer_id - 1];
    for (auto &x: r_v[layer_id]) x = F::random();
    auto previousRandom = F_ZERO;

    verify_timer.stop();
    verify_slow_timer.stop();

    p->sumcheckInitPhase2();
    for (int j = 0; j < C.circuit[layer_id].maxDadBitLength; ++j)
    {
        verify_timer.start();
        verify_slow_timer.start();
        auto poly = p->sumcheckUpdatePhase2(previousRandom);
        if (poly.eval(0) + poly.eval(1) != previousSum)
        {
            cerr << (poly.eval(0) + poly.eval(1)).toString() << ' ' << previousSum.toString() << endl;
            fprintf(stderr, "Verification fail, phase2, circuit level %d, current bit %d, total is %d\n", layer_id, j, C.circuit[layer_id].maxDadBitLength);
            return false;
        }
        previousRandom = r_v[layer_id][j];
        previousSum = poly.eval(previousRandom);

        verify_timer.stop();
        verify_slow_timer.stop();
    }

    p->sumcheckFinalize2(previousRandom, final_claims_v[layer_id].begin());

    verify_slow_timer.start();

    betaInitPhase2(layer_id);
    predicatePhase2(layer_id);

    verify_slow_timer.stop();
    return true;
}

bool verifier::verifyLiu(int layer_id, F &previousSum) {
    verify_timer.start();
    verify_slow_timer.start();
    int pre_layer_id = layer_id - 1;
    auto &cur = C.circuit[layer_id], &pre = C.circuit[pre_layer_id];

    for (auto &x: sig) x = F::random();
    for (auto &x: r_liu) x = F::random();

    previousSum = sig[0] * final_claim_u;
    for (int j = layer_id; j < C.size; ++j)
        if (~C.circuit[j].dadBitLength[pre_layer_id])
            previousSum += sig[j - pre_layer_id] * final_claims_v[j][pre_layer_id];

    verify_timer.stop();
    verify_slow_timer.stop();

    p->sumcheckInitLiu(sig.begin());
    F previousRandom = F_ZERO;
    for (u64 j = 0; j < pre.bitLength; ++j)
    {
        auto poly = p->sumcheckLiuUpdate(previousRandom);
        verify_timer.start();
        verify_slow_timer.start();
        if (poly.eval(0) + poly.eval(1) != previousSum)
        {
            cerr << (poly.eval(0) + poly.eval(1)).toString() << ' ' << previousSum.toString() << endl;
            fprintf(stderr, "Liu fail, circuit %d, current bit %d\n", layer_id, j);
            return false;
        }
        previousRandom = r_liu[j];
        previousSum = poly.eval(previousRandom);
        verify_timer.stop();
        verify_slow_timer.stop();
    }
    F gr = F_ZERO, vr;
    p->sumcheckLiuFinalize(previousRandom, vr);

    verify_slow_timer.start();
    initBetaTable(beta_u, pre.bitLength, r_liu.begin(), F_ONE);

    initBetaTable(beta_g, pre.bitLength, r_u.begin(), sig[0]);

    for (u64 g = 0; g < pre.size; ++g)
        gr = gr + beta_g[g] * beta_u[g];

    for (int j = layer_id; j < C.size; ++j)
        if (~C.circuit[j].dadBitLength[pre_layer_id]) {
            initBetaTable(beta_g, C.circuit[j].dadBitLength[pre_layer_id], r_v[j].begin(), sig[j - pre_layer_id]);
            for (u64 g = 0; g < C.circuit[j].dadSize[pre_layer_id]; ++g)
                gr = gr + beta_g[g] * beta_u[C.circuit[j].dadId[pre_layer_id][g]];
        }

    verify_timer.start();
    if (vr * gr != previousSum)
    {
        cerr << (vr * gr).toString() << ' ' << previousSum.toString() << endl;
        fprintf(stderr, "Liu fail, semi final, circuit %d\n", layer_id);
        return false;
    }
    previousSum = vr;

    verify_timer.stop();
    verify_slow_timer.stop();
    return true;
}
double verifier::verifyTime() const {
    return verify_timer.elapse_sec();
}

double verifier::verifySlowTime() const {
    return verify_slow_timer.elapse_sec();
}


#ifdef USE_VIRGO
void public_array_prepare_generic(vector<F> &q_coef_arr, vector<F> &public_array, int log_length)
{
    using namespace virgo;
	q_coef_arr.resize(1ULL << log_length);
	int coef_slice_size = (1 << (log_length - log_slice_number));
	for (int i = 0; i < (1 << log_slice_number); ++i)
	{
		inverse_fast_fourier_transform((public_array.begin() + i * coef_slice_size).base(),
                                       coef_slice_size,
                                       coef_slice_size,
                                       F::getRootOfUnity(log_length - log_slice_number),
                                       (q_coef_arr.begin() + i * coef_slice_size).base());
	}
}

bool verifier::verifyPoly(const virgo::__hhash_digest &merkle_root_l, const F &previousSum) {
    using namespace virgo;
    verify_timer.start();
    verify_slow_timer.start();

    vector<F> output(1ULL << C.circuit[0].bitLength);
    initBetaTable(output, C.circuit[0].bitLength, r_liu.begin(), F_ONE);
    vector<F> processed;
    public_array_prepare_generic(processed, output, C.circuit[0].bitLength);

    verify_timer.stop();
    verify_slow_timer.stop();

    F input_0;
    auto mask = std::vector<F>(1, F_ZERO);
    vector<F> all_sum(slice_number + 1);
    auto merkle_root_h = p->commit_public(output, input_0, mask, all_sum);
    bool flag = poly_ver.verify_poly_commitment(all_sum.data(), C.circuit[0].bitLength, processed.data(), mask, commit_vt, commit_ps, commit_pt, merkle_root_l, merkle_root_h);
    commit_ps += sizeof(__hhash_digest) * 2 + sizeof(F);

    if (previousSum != input_0 || !flag)
    {
        fprintf(stderr, "Verification fail, final input check fail.\n");
        return false;
    }
    return true;
}
#endif
