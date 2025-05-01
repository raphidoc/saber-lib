//
// Created by raphael on 5/1/25.
//

#ifndef SABER_LIB_SABER_DATA_H
#define SABER_LIB_SABER_DATA_H

#include <stddef.h>

// Global memory setters
int load_pure_water(const double* wl, const double* a, size_t n);
int load_a0_a1(const double* wl, const double* a0, const double* a1, size_t n);
int load_bottom_reflectance(const double* wl, const double* matrix, const char** colnames, size_t n_wl, size_t n_classes);

// Cache builder
int build_cache(const double* wl, size_t n);

// Cached data getters
const double* get_a_w();
const double* get_bb_w();
const double* get_a0();
const double* get_a1();
const double* get_r_rs_b();
const char**  get_r_rs_b_colnames();
size_t get_n();
size_t get_n_classes();

#endif //SABER_LIB_SABER_DATA_H
