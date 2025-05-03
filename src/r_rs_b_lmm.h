//
// Created by raphael on 5/1/25.
//

#ifndef SABER_LIB_R_RS_B_LMM_H
#define SABER_LIB_R_RS_B_LMM_H

#ifndef SABER_RRS_H
#define SABER_RRS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int compute_r_rs_b_lmm(
        const char** class_names, const double* class_fractions, size_t n_frac,
        double* out_r_rs_b  // output: array of length saber_get_n()
);

#ifdef __cplusplus
}
#endif

#endif


#endif //SABER_LIB_R_RS_B_LMM_H
