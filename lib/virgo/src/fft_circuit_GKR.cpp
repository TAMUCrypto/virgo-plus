#include <cstdio>

#include "fieldElement.hpp"
#include <vector>
#include "polynomial.h"
#include <chrono>
namespace virgo {
    namespace fft_circuit_gkr {
        class circuit {
        public:
            std::vector<fieldElement *> circuit_val;
            std::vector<int> size;
        } C;


        fieldElement rou, inv_rou;
        fieldElement inv_n;
        fieldElement *eval_points;
        int proof_size = 0;
        float v_time = 0;

        void build_circuit(int lg_size, fieldElement *r) {
            C.circuit_val.push_back(new fieldElement[1]);
            C.circuit_val[0][0] = fieldElement(1); //input layer
            C.size.push_back(1);
            for (int i = 0; i < lg_size; ++i) {
                C.circuit_val.push_back(new fieldElement[1 << (i + 1)]);
                for (int j = 0; j < (1 << i); ++j) {
                    C.circuit_val[i + 1][j << 1] = C.circuit_val[i][j] * r[i];
                    C.circuit_val[i + 1][j << 1 | 1] = C.circuit_val[i][j] * (fieldElement(1) - r[i]);
                }
                C.size.push_back(1 << (i + 1));
            }
            //ifft
            rou = fieldElement::getRootOfUnity(lg_size);
            inv_rou = rou.inv();
            inv_n = fieldElement::fastPow(fieldElement((1 << lg_size)), fieldElement::mod - 2);
            fieldElement *x_arr = new fieldElement[1 << lg_size];
            fieldElement *rot_mul = new fieldElement[fieldElement::__max_order];
            rot_mul[0] = inv_rou;
            for (int i = 1; i < fieldElement::__max_order; ++i) {
                rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1];
            }
            int starting_depth = lg_size;
            for (int dep = lg_size - 1; dep >= 0; --dep) {
                int blk_size = 1 << (lg_size - dep);
                int half_blk_size = blk_size >> 1;
                int cur = starting_depth + (lg_size - dep);
                int pre = cur - 1;
                C.circuit_val.push_back(new fieldElement[1 << lg_size]);
                C.size.push_back(1 << lg_size);
                fieldElement x = fieldElement(1);
                {
                    x_arr[0] = fieldElement(1);
                    for (int j = 1; j < blk_size; ++j)
                        x_arr[j] = x_arr[j - 1] * rot_mul[dep];
                    for (int k = 0; k < blk_size / 2; ++k) {
                        int double_k = (k) & (half_blk_size - 1);
                        for (int j = 0; j < (1 << dep); ++j) {
                            auto l_value = C.circuit_val[pre][double_k << (dep + 1) | j], r_value =
                                    x_arr[k] * C.circuit_val[pre][double_k << (dep + 1) | (1 << dep) | j];
                            C.circuit_val[cur][k << dep | j] = l_value + r_value;
                            C.circuit_val[cur][(k + blk_size / 2) << dep | j] = l_value - r_value;
                        }
                    }
                }
            }
            int last_layer_id = C.circuit_val.size();
            C.circuit_val.push_back(new fieldElement[1 << lg_size]);
            C.size.push_back(1 << lg_size);
            for (int i = 0; i < (1 << lg_size); ++i) {
                C.circuit_val[last_layer_id][i] = C.circuit_val[last_layer_id - 1][i] * inv_n;
            }

            //poly eval assuming 64 repetitions
            //directly evaluate the polynomial for practical purpose
            //since the number of repetition is only 40
            fieldElement eval_point;
            last_layer_id = C.circuit_val.size();
            C.circuit_val.push_back(new fieldElement[64 * (1 << lg_size)]);
            C.size.push_back(64 << lg_size);
            eval_points = new fieldElement[64];
            for (int i = 0; i < 64; ++i) {
                eval_point = fieldElement::random();
                eval_points[i] = eval_point;
                fieldElement x = fieldElement(1);
                for (int j = 0; j < (1 << lg_size); ++j) {
                    C.circuit_val[last_layer_id][j + (i << lg_size)] = C.circuit_val[last_layer_id - 1][j] * x;
                    x = x * eval_point;
                }
            }
            last_layer_id = C.circuit_val.size();
            C.circuit_val.push_back(new fieldElement[64]);
            C.size.push_back(64);
            for (int i = 0; i < 64; ++i) {
                C.circuit_val[last_layer_id][i] = fieldElement(0);
                for (int j = 0; j < (1 << lg_size); ++j) {
                    C.circuit_val[last_layer_id][i] =
                            C.circuit_val[last_layer_id][i] + C.circuit_val[last_layer_id - 1][j + i * (1 << lg_size)];
                }
            }
        }

        void refresh_randomness(fieldElement *r, fieldElement *one_minus_r, int size) {
            for (int i = 0; i < size; ++i) {
                r[i] = fieldElement::random();
                one_minus_r[i] = fieldElement(1) - r[i];
            }
        }

        int mylog(long long x) {
            for (int i = 0; i < 64; ++i) {
                if ((1LL << i) == x)
                    return i;
            }
            assert(false);
            return 0;
        }

        fieldElement
        V_output(const fieldElement *one_minus_r_0, const fieldElement *r_0, const fieldElement *output_raw,
                 int r_0_size, int output_size) {
            fieldElement *output;
            output = new fieldElement[output_size];
            for (int i = 0; i < output_size; ++i)
                output[i] = output_raw[i];
            for (int i = 0; i < r_0_size; ++i) {
                for (int j = 0; j < (output_size >> 1); ++j) {
                    output[j] = (output[j << 1] * one_minus_r_0[i] + output[j << 1 | 1] * r_0[i]);
                }
                output_size >>= 1;
            }
            fieldElement res = output[0];
            delete[] output;
            return res;
        }

        fieldElement alpha, beta;
        fieldElement *beta_g_r0_fhalf, *beta_g_r0_shalf, *beta_g_r1_fhalf, *beta_g_r1_shalf, *beta_u_fhalf, *beta_u_shalf;
        linear_poly *add_mult_sum, *V_mult_add, *addV_array;

        void init_array(int max_bit_length) {
            add_mult_sum = new linear_poly[(1 << max_bit_length)];
            V_mult_add = new linear_poly[(1 << max_bit_length)];
            addV_array = new linear_poly[(1 << max_bit_length)];
            int half_length = (max_bit_length >> 1) + 1;
            beta_g_r0_fhalf = new fieldElement[(1 << half_length)];
            beta_g_r0_shalf = new fieldElement[(1 << half_length)];
            beta_g_r1_fhalf = new fieldElement[(1 << half_length)];
            beta_g_r1_shalf = new fieldElement[(1 << half_length)];
            beta_u_fhalf = new fieldElement[(1 << half_length)];
            beta_u_shalf = new fieldElement[(1 << half_length)];
        }

        quadratic_poly sumcheck_phase1_update(fieldElement previous_random, int current_bit, int total_uv) {
            quadratic_poly ret = quadratic_poly(fieldElement(0), fieldElement(0), fieldElement(0));
            for (int i = 0; i < (total_uv >> 1); ++i) {
                fieldElement zero_value, one_value;
                int g_zero = i << 1, g_one = i << 1 | 1;
                if (current_bit == 0) {
                    V_mult_add[i].b = V_mult_add[g_zero].b;
                    V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

                    addV_array[i].b = addV_array[g_zero].b;
                    addV_array[i].a = addV_array[g_one].b - addV_array[i].b;

                    add_mult_sum[i].b = add_mult_sum[g_zero].b;
                    add_mult_sum[i].a = add_mult_sum[g_one].b - add_mult_sum[i].b;

                } else {
                    V_mult_add[i].b = (V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b);
                    V_mult_add[i].a = (V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b - V_mult_add[i].b);

                    addV_array[i].b = (addV_array[g_zero].a * previous_random + addV_array[g_zero].b);
                    addV_array[i].a = (addV_array[g_one].a * previous_random + addV_array[g_one].b - addV_array[i].b);

                    add_mult_sum[i].b = (add_mult_sum[g_zero].a * previous_random + add_mult_sum[g_zero].b);
                    add_mult_sum[i].a = (add_mult_sum[g_one].a * previous_random + add_mult_sum[g_one].b -
                                         add_mult_sum[i].b);

                }
                ret.a = (ret.a + add_mult_sum[i].a * V_mult_add[i].a);
                ret.b = (ret.b + add_mult_sum[i].a * V_mult_add[i].b + add_mult_sum[i].b * V_mult_add[i].a
                         + addV_array[i].a);
                ret.c = (ret.c + add_mult_sum[i].b * V_mult_add[i].b
                         + addV_array[i].b);
            }
            return ret;
        }

        quadratic_poly sumcheck_phase2_update(fieldElement previous_random, int current_bit, int total_uv) {
            quadratic_poly ret = quadratic_poly(fieldElement(0), fieldElement(0), fieldElement(0));
            for (int i = 0; i < (total_uv >> 1); ++i) {
                int g_zero = i << 1, g_one = i << 1 | 1;
                if (current_bit == 0) {
                    V_mult_add[i].b = V_mult_add[g_zero].b;
                    V_mult_add[i].a = V_mult_add[g_one].b - V_mult_add[i].b;

                    addV_array[i].b = addV_array[g_zero].b;
                    addV_array[i].a = addV_array[g_one].b - addV_array[i].b;

                    add_mult_sum[i].b = add_mult_sum[g_zero].b;
                    add_mult_sum[i].a = add_mult_sum[g_one].b - add_mult_sum[i].b;
                } else {

                    V_mult_add[i].b = (V_mult_add[g_zero].a * previous_random + V_mult_add[g_zero].b);
                    V_mult_add[i].a = (V_mult_add[g_one].a * previous_random + V_mult_add[g_one].b - V_mult_add[i].b);

                    addV_array[i].b = (addV_array[g_zero].a * previous_random + addV_array[g_zero].b);
                    addV_array[i].a = (addV_array[g_one].a * previous_random + addV_array[g_one].b - addV_array[i].b);

                    add_mult_sum[i].b = (add_mult_sum[g_zero].a * previous_random + add_mult_sum[g_zero].b);
                    add_mult_sum[i].a = (add_mult_sum[g_one].a * previous_random + add_mult_sum[g_one].b -
                                         add_mult_sum[i].b);
                }

                ret.a = (ret.a + add_mult_sum[i].a * V_mult_add[i].a);
                ret.b = (ret.b + add_mult_sum[i].a * V_mult_add[i].b
                         + add_mult_sum[i].b * V_mult_add[i].a
                         + addV_array[i].a);
                ret.c = (ret.c + add_mult_sum[i].b * V_mult_add[i].b
                         + addV_array[i].b);
            }
            return ret;
        }

        bool addition_layer(fieldElement &alpha_beta_sum, fieldElement *c_val, int size_poly, int num_poly,
                            fieldElement *r_0, fieldElement *one_minus_r_0, fieldElement *r_1,
                            fieldElement *one_minus_r_1, fieldElement *r_u, fieldElement *one_minus_r_u,
                            fieldElement *r_v, fieldElement *one_minus_r_v) {
            /*
             * addition layer init phase 1
             */
            fieldElement zero = fieldElement(0);
            beta_g_r0_fhalf[0] = alpha;
            beta_g_r1_fhalf[0] = beta;
            beta_g_r0_shalf[0] = fieldElement(1);
            beta_g_r1_shalf[0] = fieldElement(1);
            int length_g = mylog(64);
            int total_uv = size_poly * num_poly;
            int first_half = length_g >> 1, second_half = length_g - first_half;

            for (int i = 0; i < first_half; ++i) {
                for (int j = 0; j < (1 << i); ++j) {
                    beta_g_r0_fhalf[j | (1 << i)] = beta_g_r0_fhalf[j] * r_0[i];
                    beta_g_r0_fhalf[j] = beta_g_r0_fhalf[j] * one_minus_r_0[i];
                    beta_g_r1_fhalf[j | (1 << i)] = beta_g_r1_fhalf[j] * r_1[i];
                    beta_g_r1_fhalf[j] = beta_g_r1_fhalf[j] * one_minus_r_1[i];
                }
            }

            for (int i = 0; i < second_half; ++i) {
                for (int j = 0; j < (1 << i); ++j) {
                    beta_g_r0_shalf[j | (1 << i)] = beta_g_r0_shalf[j] * r_0[i + first_half];
                    beta_g_r0_shalf[j] = beta_g_r0_shalf[j] * one_minus_r_0[i + first_half];
                    beta_g_r1_shalf[j | (1 << i)] = beta_g_r1_shalf[j] * r_1[i + first_half];
                    beta_g_r1_shalf[j] = beta_g_r1_shalf[j] * one_minus_r_1[i + first_half];
                }
            }
            for (int i = 0; i < total_uv; ++i) {
                V_mult_add[i] = c_val[i];

                addV_array[i].a = zero;
                addV_array[i].b = zero;
                add_mult_sum[i].a = zero;
                add_mult_sum[i].b = zero;
            }

            int mask_fhalf = (1 << first_half) - 1;
            for (int i = 0; i < num_poly; ++i) {
                auto tmp = (beta_g_r0_fhalf[i & mask_fhalf] * beta_g_r0_shalf[i >> first_half]
                            + beta_g_r1_fhalf[i & mask_fhalf] * beta_g_r1_shalf[i >> first_half]);
                for (int j = i * size_poly; j < (i + 1) * size_poly; ++j) {
                    add_mult_sum[j].b = (add_mult_sum[j].b + tmp);
                }
            }

            //update phase 1
            int log_uv = mylog(size_poly) + mylog(num_poly);
            fieldElement previous_random = fieldElement(0);
            refresh_randomness(r_u, one_minus_r_u, log_uv);
            refresh_randomness(r_v, one_minus_r_v, log_uv);
            for (int i = 0; i < log_uv; ++i) {
                quadratic_poly poly = sumcheck_phase1_update(previous_random, i, (1 << (log_uv - i)));
                proof_size += sizeof(quadratic_poly);
                previous_random = r_u[i];
                if (poly.eval(0) + poly.eval(1) != alpha_beta_sum) {
                    return 0;
                }
                alpha_beta_sum = poly.eval(r_u[i]);
            }

            auto v_u = V_mult_add[0].eval(previous_random);
            first_half = log_uv >> 1, second_half = log_uv - first_half;

            mask_fhalf = (1 << first_half) - 1;
            int first_g_half = (length_g >> 1);
            int mask_g_fhalf = (1 << (length_g >> 1)) - 1;

            //Verifier's computation

            std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            fieldElement summation_val = fieldElement(0);
            int log_g = mylog(num_poly);
            for (int i = 0; i < num_poly; ++i) {
                auto tmp_g = (beta_g_r0_fhalf[i & mask_g_fhalf] * beta_g_r0_shalf[i >> first_g_half]
                              + beta_g_r1_fhalf[i & mask_g_fhalf] * beta_g_r1_shalf[i >> first_g_half]);
                auto tmp_u = fieldElement(1);
                for (int j = 0; j < log_g; ++j) {
                    if ((i >> j) & 1) {
                        tmp_u = tmp_u * r_u[log_uv - log_g + j];
                    } else
                        tmp_u = tmp_u * (fieldElement(1) - r_u[log_uv - log_g + j]);
                }
                summation_val = summation_val + tmp_g * tmp_u;
            }
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
                    t1 - t0);
            v_time += time_span.count();

            for (int i = 0; i < log_uv; ++i) {
                r_0[i] = r_u[i];
                r_1[i] = r_v[i];
                one_minus_r_0[i] = one_minus_r_u[i];
                one_minus_r_1[i] = one_minus_r_v[i];
            }
            if (alpha_beta_sum != summation_val * v_u)
                return 0;
            alpha_beta_sum = alpha * v_u;
            return 1;
        }

        bool
        mult_layer(fieldElement &alpha_beta_sum, fieldElement *c_val, int size_poly, int num_poly, fieldElement *r_0,
                   fieldElement *one_minus_r_0, fieldElement *r_1, fieldElement *one_minus_r_1, fieldElement *r_u,
                   fieldElement *one_minus_r_u, fieldElement *r_v, fieldElement *one_minus_r_v) {
            /*
             * mult layer init phase 1
             */
            fieldElement zero = fieldElement(0);
            beta_g_r0_fhalf[0] = alpha;
            beta_g_r1_fhalf[0] = beta;
            beta_g_r0_shalf[0] = fieldElement(1);
            beta_g_r1_shalf[0] = fieldElement(1);
            int length_g = mylog(size_poly * num_poly);
            int total_uv = size_poly;
            int first_half = length_g >> 1, second_half = length_g - first_half;

            for (int i = 0; i < first_half; ++i) {
                for (int j = 0; j < (1 << i); ++j) {
                    beta_g_r0_fhalf[j | (1 << i)] = beta_g_r0_fhalf[j] * r_0[i];
                    beta_g_r0_fhalf[j] = beta_g_r0_fhalf[j] * one_minus_r_0[i];
                    beta_g_r1_fhalf[j | (1 << i)] = beta_g_r1_fhalf[j] * r_1[i];
                    beta_g_r1_fhalf[j] = beta_g_r1_fhalf[j] * one_minus_r_1[i];
                }
            }

            for (int i = 0; i < second_half; ++i) {
                for (int j = 0; j < (1 << i); ++j) {
                    beta_g_r0_shalf[j | (1 << i)] = beta_g_r0_shalf[j] * r_0[i + first_half];
                    beta_g_r0_shalf[j] = beta_g_r0_shalf[j] * one_minus_r_0[i + first_half];
                    beta_g_r1_shalf[j | (1 << i)] = beta_g_r1_shalf[j] * r_1[i + first_half];
                    beta_g_r1_shalf[j] = beta_g_r1_shalf[j] * one_minus_r_1[i + first_half];
                }
            }
            for (int i = 0; i < total_uv; ++i) {
                V_mult_add[i] = c_val[i];

                addV_array[i].a = zero;
                addV_array[i].b = zero;
                add_mult_sum[i].a = zero;
                add_mult_sum[i].b = zero;
            }

            int mask_fhalf = (1 << first_half) - 1;
            fieldElement *x_eval = new fieldElement[num_poly];
            for (int i = 0; i < num_poly; ++i)
                x_eval[i] = fieldElement(1);
            for (int i = 0; i < size_poly; ++i) {
                for (int j = 0; j < num_poly; ++j) {
                    int g_id = j * size_poly + i;
                    auto tmp = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                                + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
                    add_mult_sum[i].b = (add_mult_sum[i].b + tmp * x_eval[j]);
                    x_eval[j] = x_eval[j] * eval_points[j];
                }
            }

            //update phase 1
            int log_uv = mylog(size_poly);
            fieldElement previous_random = fieldElement(0);
            refresh_randomness(r_u, one_minus_r_u, log_uv);
            refresh_randomness(r_v, one_minus_r_v, log_uv);
            for (int i = 0; i < log_uv; ++i) {
                quadratic_poly poly = sumcheck_phase1_update(previous_random, i, (1 << (log_uv - i)));
                proof_size += sizeof(quadratic_poly);
                previous_random = r_u[i];
                if (poly.eval(0) + poly.eval(1) != alpha_beta_sum) {
                    return 0;
                }
                alpha_beta_sum = poly.eval(r_u[i]);
            }
            auto v_u = V_mult_add[0].eval(previous_random);

            //Verifier's computation

            std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            auto summation_mult = fieldElement(0);
            int log_g = length_g;
            for (int i = 0; i < num_poly; ++i) {
                auto tmp_g_0 = alpha;
                auto tmp_g_1 = beta;
                int lg_num_poly = mylog(num_poly);
                for (int j = 0; j < lg_num_poly; ++j) {
                    if ((i >> j) & 1) {
                        tmp_g_0 = tmp_g_0 * r_0[log_g - lg_num_poly + j];
                        tmp_g_1 = tmp_g_1 * r_1[log_g - lg_num_poly + j];
                    } else {
                        tmp_g_0 = tmp_g_0 * (fieldElement(1) - r_0[log_g - lg_num_poly + j]);
                        tmp_g_1 = tmp_g_1 * (fieldElement(1) - r_1[log_g - lg_num_poly + j]);
                    }
                }
                auto tmp_u_0 = fieldElement(1), tmp_u_1 = fieldElement(1);
                auto x = eval_points[i];
                for (int j = 0; j < log_uv; ++j) {
                    tmp_u_0 = tmp_u_0 * (r_0[j] * r_u[j] * x + (fieldElement(1) - r_0[j]) * (fieldElement(1) - r_u[j]));
                    tmp_u_1 = tmp_u_1 * (r_1[j] * r_u[j] * x + (fieldElement(1) - r_1[j]) * (fieldElement(1) - r_u[j]));
                    x = x * x;
                }
                summation_mult = summation_mult + tmp_g_0 * tmp_u_0 + tmp_g_1 * tmp_u_1;
            }
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
                    t1 - t0);
            v_time += time_span.count();
            for (int i = 0; i < log_uv; ++i) {
                r_0[i] = r_u[i];
                r_1[i] = r_v[i];
                one_minus_r_0[i] = one_minus_r_u[i];
                one_minus_r_1[i] = one_minus_r_v[i];
            }
            if (alpha_beta_sum != summation_mult * v_u)
                return 0;
            alpha_beta_sum = alpha * v_u;
            return 1;
        }

        bool intermediate_layer(fieldElement &alpha_beta_sum, fieldElement *c_val, int size_poly, fieldElement *r_0,
                                fieldElement *one_minus_r_0, fieldElement *r_1, fieldElement *one_minus_r_1,
                                fieldElement *r_u, fieldElement *one_minus_r_u, fieldElement *r_v,
                                fieldElement *one_minus_r_v) {
            alpha_beta_sum = alpha_beta_sum * fieldElement(size_poly);
            assert(fieldElement(size_poly) *inv_n == fieldElement(1));
            return 1;
        }

        bool ifft_gkr(fieldElement &alpha_beta_sum, int size_poly, fieldElement *r_0, fieldElement *one_minus_r_0,
                      fieldElement *r_1, fieldElement *one_minus_r_1, fieldElement *r_u, fieldElement *one_minus_r_u,
                      fieldElement *r_v, fieldElement *one_minus_r_v) {
            int lg_size = mylog(size_poly);
            int starting_depth = lg_size;
            inv_rou = rou.inv();
            fieldElement *rot_mul = new fieldElement[fieldElement::__max_order];
            rot_mul[0] = inv_rou;
            for (int i = 1; i < fieldElement::__max_order; ++i) {
                rot_mul[i] = rot_mul[i - 1] * rot_mul[i - 1];
            }
            for (int dep = 0; dep < lg_size; ++dep) {
                /*
                 * Circuit Description:
                        int blk_size = 1 << (lg_size - dep);
                        int half_blk_size = blk_size >> 1;
                        int cur = starting_depth + (lg_size - dep);
                        int pre = cur - 1;
                        C.circuit_val.push_back(new fieldElement[1 << lg_size]);
                        C.size.push_back(1 << lg_size);
                        fieldElement x = fieldElement(1);
                        {
                            x_arr[0] = fieldElement(1);
                            for(int j = 1; j < blk_size; ++j)
                                x_arr[j] = x_arr[j - 1] * rot_mul[dep];
                            for(int k = 0; k < blk_size / 2; ++k)
                            {
                                int double_k = (k) & (half_blk_size - 1);
                                for(int j = 0; j < (1 << dep); ++j)
                                {
                                    auto l_value = C.circuit_val[pre][double_k << (dep + 1) | j], r_value = x_arr[k] * C.circuit_val[pre][double_k << (dep + 1) | (1 << dep) | j];
                                    C.circuit_val[cur][k << dep | j] = l_value + r_value;
                                    C.circuit_val[cur][(k + blk_size / 2) << dep | j] = l_value - r_value;
                                }
                            }
                        }
                 */
                int blk_size = 1 << (lg_size - dep);
                int half_blk_size = blk_size >> 1;
                int cur = starting_depth + (lg_size - dep);
                int pre = cur - 1;
                {
                    fieldElement zero = fieldElement(0);
                    beta_g_r0_fhalf[0] = alpha;
                    beta_g_r1_fhalf[0] = beta;
                    beta_g_r0_shalf[0] = fieldElement(1);
                    beta_g_r1_shalf[0] = fieldElement(1);
                    int length_g = mylog(size_poly);
                    int total_uv = size_poly;
                    int first_half = length_g >> 1, second_half = length_g - first_half;

                    for (int i = 0; i < first_half; ++i) {
                        for (int j = 0; j < (1 << i); ++j) {
                            beta_g_r0_fhalf[j | (1 << i)] = beta_g_r0_fhalf[j] * r_0[i];
                            beta_g_r0_fhalf[j] = beta_g_r0_fhalf[j] * one_minus_r_0[i];
                            beta_g_r1_fhalf[j | (1 << i)] = beta_g_r1_fhalf[j] * r_1[i];
                            beta_g_r1_fhalf[j] = beta_g_r1_fhalf[j] * one_minus_r_1[i];
                        }
                    }

                    for (int i = 0; i < second_half; ++i) {
                        for (int j = 0; j < (1 << i); ++j) {
                            beta_g_r0_shalf[j | (1 << i)] = beta_g_r0_shalf[j] * r_0[i + first_half];
                            beta_g_r0_shalf[j] = beta_g_r0_shalf[j] * one_minus_r_0[i + first_half];
                            beta_g_r1_shalf[j | (1 << i)] = beta_g_r1_shalf[j] * r_1[i + first_half];
                            beta_g_r1_shalf[j] = beta_g_r1_shalf[j] * one_minus_r_1[i + first_half];
                        }
                    }
                    for (int i = 0; i < total_uv; ++i) {
                        V_mult_add[i] = C.circuit_val[pre][i];

                        addV_array[i].a = zero;
                        addV_array[i].b = zero;
                        add_mult_sum[i].a = zero;
                        add_mult_sum[i].b = zero;
                    }

                    int mask_fhalf = (1 << first_half) - 1;
                    fieldElement x = fieldElement(1);
                    for (int k = 0; k < blk_size / 2; ++k) {
                        for (int j = 0; j < (1 << dep); ++j) {
                            int g_id = k << (dep) | j;
                            int u_id = k << (dep + 1) | j;
                            int v_id = k << (dep + 1) | (1 << dep) | j;
                            auto tmp = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                                        + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);

                            add_mult_sum[u_id].b = add_mult_sum[u_id].b + tmp;
                            addV_array[u_id].b = addV_array[u_id].b + tmp * x * C.circuit_val[pre][v_id];


                            g_id = (k + blk_size / 2) << (dep) | j;
                            u_id = k << (dep + 1) | j;
                            v_id = k << (dep + 1) | (1 << dep) | j;
                            tmp = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                                   + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
                            add_mult_sum[u_id].b = add_mult_sum[u_id].b + tmp;
                            addV_array[u_id].b = addV_array[u_id].b - tmp * x * C.circuit_val[pre][v_id];
                        }
                        x = x * rot_mul[dep];
                    }

                    //update phase 1
                    int log_uv = mylog(size_poly);
                    fieldElement previous_random = fieldElement(0);
                    refresh_randomness(r_u, one_minus_r_u, log_uv);
                    refresh_randomness(r_v, one_minus_r_v, log_uv);
                    for (int i = 0; i < log_uv; ++i) {
                        quadratic_poly poly = sumcheck_phase1_update(previous_random, i, (1 << (log_uv - i)));
                        proof_size += sizeof(quadratic_poly);
                        previous_random = r_u[i];
                        if (poly.eval(0) + poly.eval(1) != alpha_beta_sum) {
                            return 0;
                        }
                        alpha_beta_sum = poly.eval(r_u[i]);
                    }
                    auto v_u = V_mult_add[0].eval(previous_random);

                    //phase 2 init
                    beta_u_fhalf[0] = fieldElement(1);
                    beta_u_shalf[0] = fieldElement(1);
                    for (int i = 0; i < first_half; ++i) {
                        for (int j = 0; j < (1 << i); ++j) {
                            beta_u_fhalf[j | (1 << i)] = beta_u_fhalf[j] * r_u[i];
                            beta_u_fhalf[j] = beta_u_fhalf[j] * one_minus_r_u[i];
                        }
                    }

                    for (int i = 0; i < second_half; ++i) {
                        for (int j = 0; j < (1 << i); ++j) {
                            beta_u_shalf[j | (1 << i)] = beta_u_shalf[j] * r_u[i + first_half];
                            beta_u_shalf[j] = beta_u_shalf[j] * one_minus_r_u[i + first_half];
                        }
                    }
                    for (int i = 0; i < total_uv; ++i) {
                        V_mult_add[i] = C.circuit_val[pre][i];

                        addV_array[i].a = zero;
                        addV_array[i].b = zero;
                        add_mult_sum[i].a = zero;
                        add_mult_sum[i].b = zero;
                    }

                    x = fieldElement(1);
                    for (int k = 0; k < blk_size / 2; ++k) {
                        for (int j = 0; j < (1 << dep); ++j) {
                            int g_id = k << (dep) | j;
                            int u_id = k << (dep + 1) | j;
                            int v_id = k << (dep + 1) | (1 << dep) | j;
                            auto tmp_g = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                                          + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
                            auto tmp_u = (beta_u_fhalf[u_id & mask_fhalf] * beta_u_shalf[u_id >> first_half]);
                            auto tmp_g_u = tmp_g * tmp_u;
                            add_mult_sum[v_id].b = add_mult_sum[v_id].b + tmp_g_u * x;
                            addV_array[v_id].b = addV_array[v_id].b + tmp_g_u * v_u;


                            g_id = (k + blk_size / 2) << (dep) | j;
                            u_id = k << (dep + 1) | j;
                            v_id = k << (dep + 1) | (1 << dep) | j;
                            tmp_g = (beta_g_r0_fhalf[g_id & mask_fhalf] * beta_g_r0_shalf[g_id >> first_half]
                                     + beta_g_r1_fhalf[g_id & mask_fhalf] * beta_g_r1_shalf[g_id >> first_half]);
                            tmp_u = (beta_u_fhalf[u_id & mask_fhalf] * beta_u_shalf[u_id >> first_half]);
                            tmp_g_u = tmp_g * tmp_u;
                            add_mult_sum[v_id].b = add_mult_sum[v_id].b - tmp_g_u * x;
                            addV_array[v_id].b = addV_array[v_id].b + tmp_g_u * v_u;
                        }
                        x = x * rot_mul[dep];
                    }


                    //phase 2 update, I think phase 2 can be removed
                    previous_random = fieldElement(0);
                    for (int i = 0; i < log_uv; ++i) {
                        quadratic_poly poly = sumcheck_phase2_update(previous_random, i, (1 << (log_uv - i)));
                        proof_size += sizeof(quadratic_poly);
                        previous_random = r_v[i];
                        if (poly.eval(0) + poly.eval(1) != alpha_beta_sum) {
                            return 0;
                        }
                        alpha_beta_sum = poly.eval(r_v[i]);
                    }

                    auto v_v = V_mult_add[0].eval(previous_random);


                    //verifier's computation


                    std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
                    auto summation_beta_u_0_A = fieldElement(1), summation_beta_v_0_A = fieldElement(1);
                    auto summation_beta_u_1_A = fieldElement(1), summation_beta_v_1_A = fieldElement(1);
                    x = rot_mul[dep];
                    int log_k = mylog(blk_size / 2);
                    int log_j = dep;
                    summation_beta_u_0_A =
                            (fieldElement(1) - r_0[log_uv - 1]) * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * alpha;
                    summation_beta_u_1_A =
                            (fieldElement(1) - r_1[log_uv - 1]) * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * beta;
                    summation_beta_v_0_A =
                            (fieldElement(1) - r_0[log_uv - 1]) * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * alpha;
                    summation_beta_v_1_A =
                            (fieldElement(1) - r_1[log_uv - 1]) * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * beta;

                    auto summation_beta_u_0_B = r_0[log_uv - 1] * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * alpha;
                    auto summation_beta_u_1_B = r_1[log_uv - 1] * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * beta;
                    auto summation_beta_v_0_B = r_0[log_uv - 1] * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * alpha;
                    auto summation_beta_v_1_B = r_1[log_uv - 1] * (fieldElement(1) - r_u[log_j]) * r_v[log_j] * beta;
                    for (int i = 0; i < log_k; ++i) {
                        summation_beta_u_0_A = summation_beta_u_0_A *
                                               (r_0[log_j + i] * r_u[log_j + 1 + i] * r_v[log_j + 1 + i] +
                                                (fieldElement(1) - r_0[log_j + i]) *
                                                (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                (fieldElement(1) - r_v[log_j + 1 + i]));
                        summation_beta_u_1_A = summation_beta_u_1_A *
                                               (r_1[log_j + i] * r_u[log_j + 1 + i] * r_v[log_j + 1 + i] +
                                                (fieldElement(1) - r_1[log_j + i]) *
                                                (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                (fieldElement(1) - r_v[log_j + 1 + i]));
                        summation_beta_v_0_A = summation_beta_v_0_A *
                                               (r_0[log_j + i] * r_u[log_j + 1 + i] * r_v[log_j + 1 + i] * x +
                                                (fieldElement(1) - r_0[log_j + i]) *
                                                (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                (fieldElement(1) - r_v[log_j + 1 + i]));
                        summation_beta_v_1_A = summation_beta_v_1_A *
                                               (r_1[log_j + i] * r_u[log_j + 1 + i] * r_v[log_j + 1 + i] * x +
                                                (fieldElement(1) - r_1[log_j + i]) *
                                                (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                (fieldElement(1) - r_v[log_j + 1 + i]));

                        summation_beta_u_0_B = summation_beta_u_0_B *
                                               (r_0[log_j + i] * r_u[log_j + 1 + i] * r_v[log_j + 1 + i] +
                                                (fieldElement(1) - r_0[log_j + i]) *
                                                (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                (fieldElement(1) - r_v[log_j + 1 + i]));
                        summation_beta_u_1_B = summation_beta_u_1_B *
                                               (r_1[log_j + i] * r_u[log_j + 1 + i] * r_v[log_j + 1 + i] +
                                                (fieldElement(1) - r_1[log_j + i]) *
                                                (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                (fieldElement(1) - r_v[log_j + 1 + i]));
                        summation_beta_v_0_B = summation_beta_v_0_B * ((fieldElement(1) - r_0[log_j + i]) *
                                                                       (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                                       (fieldElement(1) - r_v[log_j + 1 + i]) +
                                                                       r_0[log_j + i] * r_u[log_j + 1 + i] *
                                                                       r_v[log_j + 1 + i] * x);
                        summation_beta_v_1_B = summation_beta_v_1_B * ((fieldElement(1) - r_1[log_j + i]) *
                                                                       (fieldElement(1) - r_u[log_j + 1 + i]) *
                                                                       (fieldElement(1) - r_v[log_j + 1 + i]) +
                                                                       r_1[log_j + i] * r_u[log_j + 1 + i] *
                                                                       r_v[log_j + 1 + i] * x);
                        x = x * x;
                    }
                    for (int i = 0; i < log_j; ++i) {
                        summation_beta_u_0_A = summation_beta_u_0_A * (r_0[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_0[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));
                        summation_beta_u_1_A = summation_beta_u_1_A * (r_1[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_1[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));
                        summation_beta_v_0_A = summation_beta_v_0_A * (r_0[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_0[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));
                        summation_beta_v_1_A = summation_beta_v_1_A * (r_1[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_1[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));

                        summation_beta_u_0_B = summation_beta_u_0_B * (r_0[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_0[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));
                        summation_beta_u_1_B = summation_beta_u_1_B * (r_1[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_1[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));
                        summation_beta_v_0_B = summation_beta_v_0_B * (r_0[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_0[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));
                        summation_beta_v_1_B = summation_beta_v_1_B * (r_1[i] * r_u[i] * r_v[i] +
                                                                       (fieldElement(1) - r_1[i]) *
                                                                       (fieldElement(1) - r_u[i]) *
                                                                       (fieldElement(1) - r_v[i]));
                    }

                    if (alpha_beta_sum !=
                        (summation_beta_u_0_A + summation_beta_u_1_A + summation_beta_u_0_B + summation_beta_u_1_B) *
                        v_u +
                        (summation_beta_v_0_A + summation_beta_v_1_A - summation_beta_v_0_B - summation_beta_v_1_B) *
                        v_v)
                        return 0;

                    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
                            t1 - t0);
                    v_time += time_span.count();
                    for (int i = 0; i < log_uv; ++i) {
                        r_0[i] = r_u[i];
                        r_1[i] = r_v[i];
                        one_minus_r_0[i] = one_minus_r_u[i];
                        one_minus_r_1[i] = one_minus_r_v[i];
                    }
                    alpha = fieldElement::random();
                    beta = fieldElement::random();
                    alpha_beta_sum = alpha * v_u + beta * v_v;
                }
            }
            return 1;
        }

        bool extension_gkr(fieldElement &alpha_beta_sum, int size_poly, fieldElement *r_0, fieldElement *one_minus_r_0,
                           fieldElement *r_1, fieldElement *one_minus_r_1, fieldElement *r_u,
                           fieldElement *one_minus_r_u, fieldElement *r_v, fieldElement *one_minus_r_v) {
            int lg_size = mylog(size_poly);
            for (int i = 1; i <= lg_size; ++i) {
                for (int j = 0; j < i; ++j)
                    proof_size += sizeof(quadratic_poly);
            }
            return 1;
        }

        float p_time = 0;

        int engage_gkr(int lg_size) {
            alpha = fieldElement(1);
            beta = fieldElement(0);

            //layer output, summation over 64 polynomials
            fieldElement *r_0, *one_minus_r_0, *r_1, *one_minus_r_1;
            fieldElement *r_u, *one_minus_r_u, *r_v, *one_minus_r_v;
            r_0 = new fieldElement[lg_size + 10];
            one_minus_r_0 = new fieldElement[lg_size + 10];
            r_1 = new fieldElement[lg_size + 10];
            one_minus_r_1 = new fieldElement[lg_size + 10];
            r_u = new fieldElement[lg_size + 10];
            one_minus_r_u = new fieldElement[lg_size + 10];
            r_v = new fieldElement[lg_size + 10];
            one_minus_r_v = new fieldElement[lg_size + 10];
            refresh_randomness(r_0, one_minus_r_0, lg_size + 10);
            refresh_randomness(r_1, one_minus_r_1, lg_size + 10);

            fieldElement a_0 = V_output(one_minus_r_0, r_0, C.circuit_val[C.circuit_val.size() - 1],
                                        mylog(C.size[C.size.size() - 1]), C.size[C.size.size() - 1]);
            fieldElement alpha_beta_sum = a_0;

            bool result = true;

            //poly eval part
            //last layer, only addition
            std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
            result &= addition_layer(alpha_beta_sum, C.circuit_val[C.circuit_val.size() - 2],
                                     C.size[C.size.size() - 2] / 64, 64, r_0, one_minus_r_0, r_1, one_minus_r_1, r_u,
                                     one_minus_r_u, r_v, one_minus_r_v);
            result &= mult_layer(alpha_beta_sum, C.circuit_val[C.circuit_val.size() - 3], C.size[C.size.size() - 3], 64,
                                 r_0, one_minus_r_0, r_1, one_minus_r_1, r_u, one_minus_r_u, r_v, one_minus_r_v);
            result &= intermediate_layer(alpha_beta_sum, C.circuit_val[C.circuit_val.size() - 4],
                                         C.size[C.size.size() - 4], r_0, one_minus_r_0, r_1, one_minus_r_1, r_u,
                                         one_minus_r_u, r_v, one_minus_r_v);
            //ifft part
            result &= ifft_gkr(alpha_beta_sum, C.size[C.size.size() - 4], r_0, one_minus_r_0, r_1, one_minus_r_1, r_u,
                               one_minus_r_u, r_v, one_minus_r_v);
            //extension part
            result &= extension_gkr(alpha_beta_sum, C.size[C.size.size() - 4], r_0, one_minus_r_0, r_1, one_minus_r_1,
                                    r_u, one_minus_r_u, r_v, one_minus_r_v);
            std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(
                    t1 - t0);
            p_time += time_span.count();
            return result;
        }

        int fft_gkr(int lg_size, double &vt, int &ps, double &pt) {
            v_time = 0;
            proof_size = 0;
            p_time = 0;
            init_array(lg_size + 6);
            fieldElement *r = new fieldElement[lg_size];
            for (int i = 0; i < lg_size; ++i)
                r[i] = fieldElement::random();
            build_circuit(lg_size, r); //not strict arith circuit
            bool result = engage_gkr(lg_size);
            if (!result)
                fprintf(stderr, "Error, fft gkr failed\n");
            vt = v_time;
            ps = proof_size;
            pt = p_time;
            return 0;
        }
    }
}