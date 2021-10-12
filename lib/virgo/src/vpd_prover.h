#ifndef __vpd_prover
#define __vpd_prover
#include <vector>
#include "fieldElement.hpp"
#include "my_hhash.h"
namespace virgo {
    __hhash_digest
    vpd_prover_init(fieldElement *l_eval, fieldElement *&l_coef, int log_input_length, int slice_size, int slice_count);
}
#endif