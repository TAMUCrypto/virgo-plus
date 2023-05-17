#pragma once

#include "circuit.h"
#include <vector>
#include <ctime>
#include <algorithm>
#include <memory>
#include "config_pc.hpp"
#include "polynomial.h"

using std::unique_ptr;
class prover
{
public:
	explicit prover(const layeredCircuit &cir);
	void evaluate();
	void init();
    void sumcheckInitAll(const vector<F>::const_iterator &r_last);
    void sumcheckInit();
    void sumcheckInitPhase1(const F &assert_random);
	void sumcheckInitPhase2();
    void sumcheckInitLiu(vector<F>::const_iterator s);

	quadratic_poly sumcheckUpdatePhase1(const F &previousRandom);
	quadratic_poly sumcheckUpdatePhase2(const F &previousRandom);
	quadratic_poly sumcheckLiuUpdate(const F &previousRandom);

    void sumcheckFinalize1(const F &previousRandom, F &claim);
    void sumcheckFinalize2(const F &previousRandom, vector<F>::iterator claims);
    void sumcheckLiuFinalize(const F &previousRandom, F &claim);

    F Vres(const vector<F>::const_iterator &r_0, int r_0_size);

	double proveTime() const;
	double proofSize() const;

#ifdef USE_VIRGO
	virgo::poly_commit::poly_commit_prover poly_prover;
	virgo::__hhash_digest commit_private();
	F inner_prod(const vector<F> &a, const vector<F> &b, u64 l);
	virgo::__hhash_digest commit_public(vector<F> &pub, F &inner_product_sum, std::vector<F> &mask, vector<F> &all_sum);
#endif

private:
    quadratic_poly sumcheckUpdate(const F &previous_random, vector<F> &r_arr, int n_pre_layer);
    quadratic_poly sumcheckUpdateEach(const F &previous_random, int idx);

    const layeredCircuit &C;
    vector<vector<F>> circuitValue;

    vector<F> r_u, r_liu;
    vector<vector<F>> r_v;

    vector<F> beta_g;
    vector<F> beta_u;
    vector<u64> total;
    vector<u64> totalSize;

    int round;
    int sumcheckLayerId;
    vector<vector<linear_poly>> multArray, addVArray, Vmult;
    F add_term, V_u;

    timer prove_timer;
    u64 proof_size;

};