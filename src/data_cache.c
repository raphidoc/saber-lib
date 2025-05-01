#include "../include/data_cache.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

// ---------- Global Memory ----------

static double* wl_cache = NULL;
static size_t n_cache = 0;

static double* cached_a_w = NULL;
static double* cached_a0 = NULL;
static double* cached_a1 = NULL;
static double* cached_bb_w = NULL;
static double* cached_r_rs_b = NULL;

static char**  r_rs_b_colnames = NULL;
static size_t n_classes = 0;

static double* a_w_wl = NULL, *a_w_val = NULL;
static double* a0a1_wl = NULL, *a0_val = NULL, *a1_val = NULL;
static double* r_rs_b_wl = NULL, *r_rs_b_matrix = NULL;

// ---------- Loaders ----------

int load_pure_water(const double* wl, const double* a, size_t n) {
    if (!wl || !a || n == 0) return 1;
    a_w_wl = realloc(a_w_wl, sizeof(double) * n);
    a_w_val = realloc(a_w_val, sizeof(double) * n);
    if (!a_w_wl || !a_w_val) return 2;
    memcpy(a_w_wl, wl, sizeof(double) * n);
    memcpy(a_w_val, a, sizeof(double) * n);
    return 0;
}

int load_a0_a1(const double* wl, const double* a0, const double* a1, size_t n) {
    if (!wl || !a0 || !a1 || n == 0) return 1;
    a0a1_wl = realloc(a0a1_wl, sizeof(double) * n);
    a0_val = realloc(a0_val, sizeof(double) * n);
    a1_val = realloc(a1_val, sizeof(double) * n);
    if (!a0a1_wl || !a0_val || !a1_val) return 2;
    memcpy(a0a1_wl, wl, sizeof(double) * n);
    memcpy(a0_val, a0, sizeof(double) * n);
    memcpy(a1_val, a1, sizeof(double) * n);
    return 0;
}

int load_bottom_reflectance(const double* wl, const double* matrix, const char** colnames, size_t n_wl, size_t n_cls) {
    r_rs_b_wl = realloc(r_rs_b_wl, sizeof(double) * n_wl);
    r_rs_b_matrix = realloc(r_rs_b_matrix, sizeof(double) * n_wl * n_cls);
    if (!r_rs_b_wl || !r_rs_b_matrix) return 1;

    memcpy(r_rs_b_wl, wl, sizeof(double) * n_wl);
    memcpy(r_rs_b_matrix, matrix, sizeof(double) * n_wl * n_cls);

    r_rs_b_colnames = (char**)colnames;  // assume caller persists strings
    n_classes = n_cls;
    return 0;
}

// ---------- Interpolation Functions ----------

double interpolate_scalar(const double* wl, const double* val, size_t n, double target) {
    if (target <= wl[0] || target >= wl[n-1]) return 0.0;
    for (size_t i = 0; i < n - 1; i++) {
        if (target >= wl[i] && target <= wl[i+1]) {
            double slope = (val[i+1] - val[i]) / (wl[i+1] - wl[i]);
            return val[i] + slope * (target - wl[i]);
        }
    }
    return 0.0;
}

int interpolate_vector(const double* wl_old, const double* val_old,
                       const double* wl_new, size_t n_old, size_t n_new,
                       double* val_out) {
    for (size_t i = 0; i < n_new; i++) {
        val_out[i] = interpolate_scalar(wl_old, val_old, n_old, wl_new[i]);
    }
    return 0;
}

int interpolate_matrix(const double* wl_old, const double* mat_old,
                       const double* wl_new, size_t n_old, size_t n_new, size_t n_class,
                       double* mat_out) {
    for (size_t j = 0; j < n_class; j++) {
        for (size_t i = 0; i < n_new; i++) {
            const double* col = mat_old + j * n_old;
            mat_out[i + j * n_new] = interpolate_scalar(wl_old, col, n_old, wl_new[i]);
        }
    }
    return 0;
}

// ---------- Cache Builder ----------

int build_cache(const double* wl, size_t n) {
    if (!a0a1_wl || !a0_val || !a1_val || !a_w_wl || !a_w_val || !r_rs_b_wl || !r_rs_b_matrix)
        return 1;

    wl_cache = realloc(wl_cache, sizeof(double) * n);
    memcpy(wl_cache, wl, sizeof(double) * n);
    n_cache = n;

    cached_a0 = realloc(cached_a0, sizeof(double) * n);
    cached_a1 = realloc(cached_a1, sizeof(double) * n);
    cached_a_w = realloc(cached_a_w, sizeof(double) * n);
    cached_bb_w = realloc(cached_bb_w, sizeof(double) * n);
    cached_r_rs_b = realloc(cached_r_rs_b, sizeof(double) * n * n_classes);

    interpolate_vector(a0a1_wl, a0_val, wl, n_cache, n_cache, cached_a0);
    interpolate_vector(a0a1_wl, a1_val, wl, n_cache, n_cache, cached_a1);
    interpolate_vector(a_w_wl, a_w_val, wl, n_cache, n_cache, cached_a_w);
    interpolate_matrix(r_rs_b_wl, r_rs_b_matrix, wl, n_cache, n_cache, n_classes, cached_r_rs_b);

    for (size_t i = 0; i < n; i++) {
        double lambda = wl[i];
        double b1 = 0.00111;
        double lambda1 = 500.0;
        double exponent = -4.32;
        cached_bb_w[i] = b1 * pow(lambda / lambda1, exponent);
    }

    return 0;
}

// ---------- Cache Accessors ----------

const double* get_a_w()      { return cached_a_w; }
const double* get_a0()       { return cached_a0; }
const double* get_a1()       { return cached_a1; }
const double* get_bb_w()     { return cached_bb_w; }
const double* get_r_rs_b()   { return cached_r_rs_b; }
const char**  get_r_rs_b_colnames() { return (const char**)r_rs_b_colnames; }
size_t get_n()              { return n_cache; }
size_t get_n_classes()      { return n_classes; }
