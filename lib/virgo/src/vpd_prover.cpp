#include "vpd_prover.h"
#include "RS_polynomial.h"
#include "constants.h"
#include "fri.h"

namespace virgo {
    __hhash_digest merkle_root;

    __hhash_digest vpd_prover_init(fieldElement *l_eval, fieldElement *&l_coef, int log_input_length, int slice_size,
                                   int slice_count) {
        //fft and apply mask
        merkle_root = fri::request_init_commit(log_input_length, 0);
        return merkle_root;
    }
}