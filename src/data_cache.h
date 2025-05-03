//
// Created by raphael on 5/1/25.
//

#ifndef SABER_LIB_DATA_CACHE_H
#define SABER_LIB_DATA_CACHE_H

#include <stddef.h>

// Global memory setters
int load_pure_water(const double* wl, const double* a, size_t n);
int load_a0_a1(const double* wl, const double* a0, const double* a1, size_t n);
int load_r_rs_b(const double* wl, const char** class_names, const double* matrix, size_t n_wl, size_t n_class);

// Cache builder
int build_cache(const double* wl, size_t n);
int ensure_cache(const double *wl, size_t n);
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

#endif //SABER_LIB_DATA_CACHE_H
