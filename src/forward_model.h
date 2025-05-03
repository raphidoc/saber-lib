//
// Created by raphael on 5/1/25.
//

#ifndef SABER_LIB_FORWARD_MODEL_H
#define SABER_LIB_FORWARD_MODEL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int forward_am03(
        const double* wavelength,
        const double* a,
        const double* bb,
        size_t n,
        int water_type,
        double theta_sun_deg,
        double theta_view_deg,
        int shallow,
        double h_w,
        const double* r_b,
        double* rrs_out
);

#ifdef __cplusplus
}
#endif

#endif //SABER_LIB_FORWARD_MODEL_H
