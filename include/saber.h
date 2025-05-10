#ifndef SABER_LIB_SABER_H
#define SABER_LIB_SABER_H

#include <stddef.h>
#include "saber_version.h"

//#define SABER_VERSION "0.1.2"

#ifdef __cplusplus
extern "C" {
#endif

const char *saber_version(void);

// Global memory setters
int load_pure_water(const double* wl, const double* a, size_t n);
int load_a0_a1(const double* wl, const double* a0, const double* a1, size_t n);
int load_r_rs_b(const double* wl, const char** colnames, const double* matrix, size_t wl_n, size_t class_n);

// Cache builder
int build_cache(const double* wl, size_t n);
void saber_reset_tables(void);

// Cached data getters
const double* get_a_w();
const double* get_bb_w();
const double* get_a0();
const double* get_a1();
const double* get_r_rs_b();
const char**  get_r_rs_b_class_names();
size_t get_n_wl();
size_t get_n_class();

int iop_from_oac(
        const double* wavelength, size_t n,
        const char** param_names, const double* param_values, size_t n_param,
        double* a_out, double* bb_out
);

int compute_r_rs_b_lmm(
        const char** class_names, const double* class_fractions, size_t n_frac,
        double* out_r_rs_b  // output: array of length saber_get_n()
);

void snell_law(double theta_view_deg, double theta_sun_deg,
               double* view_w, double* sun_w);

int forward_am03(
        const double *wavelength,
        const double *a,
        const double *bb,
        size_t n,
        int water_type,
        double theta_sun_deg,
        double theta_view_deg,
        int shallow,
        double h_w,
        const double *r_b,
        double *rrs_out
);

int retrieve_r_rs_b_am03(
        const double *wavelength,
        const double *a,
        const double *bb,
        const double *r_rs_obs,
        size_t n,
        int water_type,
        double theta_sun_deg,
        double theta_view_deg,
        double h_w,
        double *r_rs_b_out
);

#ifdef __cplusplus
}
#endif

#endif //SABER_LIB_SABER_H
