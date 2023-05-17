#include "RS_polynomial.h"
#include "utility.h"
#include "poly_commit.h"
#include "fft_circuit_GKR.h"
#include <iostream>
#include <cstdio>

namespace virgo {
    inline bool verify_merkle(__hhash_digest h, std::vector<__hhash_digest> merkle_path, int len, int pow,
                              std::vector<std::pair<fieldElement, fieldElement> > value) {
        __hhash_digest cur_hhash = merkle_path[len - 1];
        __hhash_digest data[2];
        for (int i = 0; i < len - 1; ++i) {
            data[0] = cur_hhash;
            data[1] = merkle_path[i];
            if ((pow & 1)) {
                data[0] = merkle_path[i];
                data[1] = cur_hhash;
            }
            pow /= 2;
            my_hhash(data, &cur_hhash);
        }
        memset(data, 0, sizeof data);
        assert(value.size() % 2 == 1);
        __hhash_digest value_h;
        memset(&value_h, 0, sizeof(__hhash_digest));
        for (int i = 0; i < value.size(); ++i) {
            fieldElement data_ele[2];
            data_ele[0] = value[i].first;
            data_ele[1] = value[i].second;
            memcpy(&data[0], data_ele, sizeof(__hhash_digest));
            data[1] = value_h;
            my_hhash(data, &value_h);
            unsigned long long *h0, *h1;
            h0 = (unsigned long long *) &value_h.h0;
            h1 = (unsigned long long *) &value_h.h1;
        }
        return equals(h, cur_hhash) && equals(value_h, merkle_path[len - 1]);
    }

    fieldElement *mask_q_coef;

//return the hhash array of commitments, randomness and final small polynomial (represented by rscode)
    poly_commit::ldt_commitment poly_commit::poly_commit_prover::commit_phase(int log_length) {
        //Prover do the work
        std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
        auto log_current_witness_size_per_slice_cp = fri::log_current_witness_size_per_slice;
        //assuming we already have the initial commit
        int codeword_size = (1 << (log_length + rs_code_rate - log_slice_number));
        //repeat until the codeword is constant
        __hhash_digest *ret = new __hhash_digest[log_length + rs_code_rate - log_slice_number];
        fieldElement *randomness = new fieldElement[log_length + rs_code_rate];
        int ptr = 0;
        while (codeword_size > (1 << rs_code_rate)) {
            assert(ptr < log_length + rs_code_rate - log_slice_number);
            randomness[ptr] = fieldElement::random();
            ret[ptr] = fri::commit_phase_step(randomness[ptr]);
            codeword_size /= 2;
            ptr++;
        }
        ldt_commitment com;
        com.commitment_hhash = ret;
        com.final_rs_code = fri::commit_phase_final();
        com.randomness = randomness;
        com.mx_depth = ptr;
        fri::log_current_witness_size_per_slice = log_current_witness_size_per_slice_cp;

        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
        total_time += time_span.count();
        //printf("Commit time %lf\n", time_span.count());

        return com;
    }

    bool poly_commit::poly_commit_verifier::verify_poly_commitment(fieldElement *all_sum, int log_length,
                                                                   fieldElement *public_array,
                                                                   std::vector<fieldElement> &all_pub_mask,
                                                                   double &v_time, int &proof_size, double &p_time,
                                                                   __hhash_digest merkle_tree_l,
                                                                   __hhash_digest merkle_tree_h) {
        auto all_pub_msk_arr = new fieldElement[all_pub_mask.size()];
        for (int i = 0; i < all_pub_mask.size(); ++i)
            all_pub_msk_arr[i] = all_pub_mask[i];
        fieldElement *all_pub_msk_coef = new fieldElement[all_pub_mask.size()];
        inverse_fast_fourier_transform(all_pub_msk_arr, all_pub_mask.size(), all_pub_mask.size(),
                                       fieldElement::getRootOfUnity(mylog(all_pub_mask.size())), all_pub_msk_coef);

        //prepare ratio and array q
        double v_time_fft, p_time_fft;
        int proof_size_fft;
        fft_circuit_gkr::fft_gkr(log_length - log_slice_number, v_time_fft, proof_size_fft, p_time_fft);
        v_time += v_time_fft;
        p_time += p_time_fft;
        proof_size += proof_size_fft;

        auto com = p->commit_phase(log_length);

        int coef_slice_size = (1 << (log_length - log_slice_number));

        for (int rep = 0; rep < 33; ++rep) {
            int slice_count = 1 << log_slice_number;
            slice_count++; //for masks
            int slice_size = (1 << (log_length + rs_code_rate - log_slice_number));

            std::chrono::high_resolution_clock::time_point t0, t1;
            std::chrono::duration<double> time_span;
            auto inv_2 = fieldElement(2).inv();
            std::pair<fieldElement, std::pair<__hhash_digest *, int> > pre_alpha_1;
            std::pair<std::vector<std::pair<fieldElement, fieldElement> >, std::vector<__hhash_digest> > alpha_l, alpha_h, alpha, beta_l, beta_h, beta;
            fieldElement s0, s1, pre_y;
            fieldElement root_of_unity;
            fieldElement y;
            bool equ_beta;
            assert(log_length - log_slice_number > 0);
            long long pow;
            for (int i = 0; i < log_length - log_slice_number; ++i) {
                t0 = std::chrono::high_resolution_clock::now();
                if (i == 0) {
                    do {
                        pow = rand() % (1 << (log_length + rs_code_rate - log_slice_number - i));
                    } while (pow < (1 << (log_length - log_slice_number - i)) || pow % 2 == 1);
                    root_of_unity = fieldElement::getRootOfUnity(log_length + rs_code_rate - log_slice_number - i);
                    y = fieldElement::fastPow(root_of_unity, pow);
                } else {
                    root_of_unity = root_of_unity * root_of_unity;
                    pow = pow % (1 << (log_length + rs_code_rate - log_slice_number - i));
                    pre_y = y;
                    y = y * y;
                }
                long long s0_pow, s1_pow;
                assert(pow % 2 == 0);
                s0_pow = pow / 2;
                s1_pow = (pow + (1LL << (log_length + rs_code_rate - log_slice_number - i))) / 2;
                s0 = fieldElement::fastPow(root_of_unity, s0_pow);
                s1 = fieldElement::fastPow(root_of_unity, s1_pow);
                int indicator;
                if (i != 0) {
                    assert(s1 == pre_y || s0 == pre_y);
                    if (s1 == pre_y)
                        indicator = 1;
                    else
                        indicator = 0;
                }
                assert(s0 * s0 == y);
                assert(s1 * s1 == y);
                int new_size = 0;
                if (i == 0) {
                    t1 = std::chrono::high_resolution_clock::now();
                    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
                    v_time += time_span.count();
                    alpha_l = fri::request_init_value_with_merkle(s0_pow, s1_pow, new_size, 0);
                    alpha_h = fri::request_init_value_with_merkle(s0_pow, s1_pow, new_size, 1);

                    proof_size += new_size; //both h and p

                    t0 = std::chrono::high_resolution_clock::now();
                    if (!verify_merkle(merkle_tree_l, alpha_l.second, alpha_l.second.size(), min(s0_pow, s1_pow),
                                       alpha_l.first))
                        return false;
                    if (!verify_merkle(merkle_tree_h, alpha_h.second, alpha_h.second.size(), min(s0_pow, s1_pow),
                                       alpha_h.first))
                        return false;
                    t1 = std::chrono::high_resolution_clock::now();
                    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
                    v_time += time_span.count();
                    beta = fri::request_step_commit(0, pow / 2, new_size);

                    proof_size += new_size;

                    t0 = std::chrono::high_resolution_clock::now();
                    if (!verify_merkle(com.commitment_hhash[0], beta.second, beta.second.size(), pow / 2, beta.first))
                        return false;

                    auto inv_mu = fieldElement::fastPow(root_of_unity, pow / 2).inv();
                    alpha.first.clear(), alpha.second.clear();
                    int mask_position_gap = slice_size / all_pub_mask.size();
                    fieldElement rou[2], x[2], inv_x[2], msk_rou[2];
                    x[0] = fieldElement::getRootOfUnity(mylog(slice_size));
                    x[1] = fieldElement::getRootOfUnity(mylog(slice_size));
                    x[0] = fieldElement::fieldElement::fastPow(x[0], s0_pow);
                    x[1] = fieldElement::fieldElement::fastPow(x[1], s1_pow);
                    msk_rou[0] = fieldElement::fieldElement::fastPow(x[0], slice_size / mask_position_gap);
                    msk_rou[1] = fieldElement::fieldElement::fastPow(x[1], slice_size / mask_position_gap);
                    rou[0] = fieldElement::fieldElement::fastPow(x[0], slice_size >> rs_code_rate);
                    rou[1] = fieldElement::fieldElement::fastPow(x[1], slice_size >> rs_code_rate);
                    inv_x[0] = x[0].inv();
                    inv_x[1] = x[1].inv();
                    fieldElement inv_msk_H;
                    fieldElement inv_H;
                    inv_msk_H = fieldElement(slice_size / mask_position_gap).inv();
                    inv_H = fieldElement(slice_size >> rs_code_rate).inv();
                    alpha.first.resize(slice_count);

                    fieldElement q_eval_0_msk, q_eval_1_msk, x0, x1, tst0, tst1;
                    fieldElement q_eval_0_val, q_eval_1_val;
                    tst0 = tst1 = q_eval_0_msk = q_eval_1_msk = fieldElement(0);
                    q_eval_0_val = q_eval_1_val = fieldElement(0);
                    x0 = x1 = fieldElement(1);
                    for (int k = 0; k < all_pub_mask.size(); ++k) {
                        q_eval_0_msk = q_eval_0_msk + x0 * all_pub_msk_coef[k];
                        x0 = x0 * x[0];
                        q_eval_1_msk = q_eval_1_msk + x1 * all_pub_msk_coef[k];
                        x1 = x1 * x[1];
                    }
                    for (int j = 0; j < slice_count; ++j) {
                        tst0 = tst1 = fieldElement(0);
                        x0 = x1 = fieldElement(1);
                        for (int k = 0; k < (1 << (log_length - log_slice_number)); ++k) {
                            if (j == slice_count - 1)
                                continue;
                            tst0 = tst0 + x0 * public_array[k + j * coef_slice_size];
                            x0 = x0 * x[0];
                            tst1 = tst1 + x1 * public_array[k + j * coef_slice_size];
                            x1 = x1 * x[1];
                        }
                        fieldElement q_eval_0 = fieldElement(0), x0 = fieldElement(1);
                        fieldElement q_eval_1 = fieldElement(0), x1 = fieldElement(1);
                        if (j != slice_count - 1) {
                            q_eval_0 = q_eval_0_val;
                            q_eval_1 = q_eval_1_val;
                        } else {
                            q_eval_0 = q_eval_0_msk;
                            q_eval_1 = q_eval_1_msk;
                        }
                        auto one = fieldElement(1);
                        //merge l and h
                        if (j == slice_count - 1) {
                            alpha.first[j].first =
                                    alpha_l.first[j].first * q_eval_0 - (msk_rou[0] - one) * alpha_h.first[j].first;
                            alpha.first[j].first =
                                    (alpha.first[j].first * fieldElement(slice_size / mask_position_gap) - all_sum[j]) *
                                    inv_x[0];
                            alpha.first[j].second =
                                    alpha_l.first[j].second * q_eval_1 - (msk_rou[1] - one) * alpha_h.first[j].second;
                            alpha.first[j].second =
                                    (alpha.first[j].second * fieldElement(slice_size / mask_position_gap) -
                                     all_sum[j]) * inv_x[1];
                        } else {
                            alpha.first[j].first =
                                    alpha_l.first[j].first * tst0 - (rou[0] - one) * alpha_h.first[j].first;
                            alpha.first[j].first =
                                    (alpha.first[j].first * fieldElement(slice_size >> rs_code_rate) - all_sum[j]) *
                                    inv_x[0];
                            alpha.first[j].second =
                                    alpha_l.first[j].second * tst1 - (rou[1] - one) * alpha_h.first[j].second;
                            alpha.first[j].second =
                                    (alpha.first[j].second * fieldElement(slice_size >> rs_code_rate) - all_sum[j]) *
                                    inv_x[1];
                        }
                        if (s0_pow > s1_pow)
                            std::swap(alpha.first[j].first, alpha.first[j].second);
                        auto p_val = (alpha.first[j].first + alpha.first[j].second) * inv_2 +
                                     (alpha.first[j].first - alpha.first[j].second) * inv_2 * com.randomness[i] *
                                     inv_mu;

                        if (p_val != beta.first[j].first && p_val != beta.first[j].second) {
                            fprintf(stderr, "Fri check consistency first round fail\n");
                            return false;
                        }
                        if (p_val == beta.first[j].first)
                            equ_beta = false;
                        else
                            equ_beta = true;
                    }

                    t1 = std::chrono::high_resolution_clock::now();
                    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
                    //This will not added into v time since the fft gkr already give the result, we didn't have time to integrate the fft gkr into the main body, so we have the evaluation code here
                    //v_time += time_span.count();
                } else {
                    t1 = std::chrono::high_resolution_clock::now();
                    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
                    v_time += time_span.count();
                    // std::cerr << "Verification Time " << v_time << std::endl;
                    if (indicator == 1) {
                        alpha = beta;
                    } else {
                        alpha = beta;
                    }

                    beta = fri::request_step_commit(i, pow / 2, new_size);

                    proof_size += new_size;

                    t0 = std::chrono::high_resolution_clock::now();
                    if (!verify_merkle(com.commitment_hhash[i], beta.second, beta.second.size(), pow / 2, beta.first))
                        return false;

                    auto inv_mu = fieldElement::fastPow(root_of_unity, pow / 2).inv();
                    t1 = std::chrono::high_resolution_clock::now();
                    time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0);
                    v_time += time_span.count();
                    for (int j = 0; j < slice_count; ++j) {
                        auto p_val_0 = (alpha.first[j].first + alpha.first[j].second) * inv_2 +
                                       (alpha.first[j].first - alpha.first[j].second) * inv_2 * com.randomness[i] *
                                       inv_mu;
                        auto p_val_1 = (alpha.first[j].first + alpha.first[j].second) * inv_2 +
                                       (alpha.first[j].second - alpha.first[j].first) * inv_2 * com.randomness[i] *
                                       inv_mu;
                        if (p_val_0 != beta.first[j].first && p_val_0 != beta.first[j].second &&
                            p_val_1 != beta.first[j].first && p_val_1 != beta.first[j].second) {
                            fprintf(stderr, "Fri check consistency %d round fail\n", i);
                            return false;
                        }
                    }
                }
            }
            delete[] pre_alpha_1.second.first;
            //CHECK last rs code
            for (int i = 0; i < slice_count - 1; ++i) {
                auto tmplate = fri::cpd.rs_codeword[com.mx_depth - 1][0 << (log_slice_number + 1) | i << 1 | 0];
                for (int j = 0; j < (1 << (rs_code_rate - 1)); ++j) {
                    if (fri::cpd.rs_codeword[com.mx_depth - 1][j << (log_slice_number + 1) | i << 1 | 0] != tmplate) {
                        fprintf(stderr, "Fri rs code check fail\n");
                        return false;
                    }
                }
            }
            for (int j = 1; j < (1 << rs_code_rate); ++j)
                if (fri::cpd.rs_codeword_msk[com.mx_depth - 1][j] !=
                    fri::cpd.rs_codeword_msk[com.mx_depth - 1][j - 1]) {
                    fprintf(stderr, "Fri msk rs code check fail\n");
                    return false;
                }
        }
        return true;
    }
}
