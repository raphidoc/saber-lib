//
// Created by raphael on 5/1/25.
//

#ifndef SABER_LIB_IOP_FROM_OAC_H
#define SABER_LIB_IOP_FROM_OAC_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int iop_from_oac(
        const double* wavelength, size_t n,
        const char** param_names, const double* param_values, size_t n_param,
        double* a_out, double* bb_out
);

#ifdef __cplusplus
}
#endif

#endif //SABER_LIB_IOP_FROM_OAC_H
