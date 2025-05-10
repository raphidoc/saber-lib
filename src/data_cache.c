#include "data_cache.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>

/* -------------------------------------------------- *
 *  64‑bit FNV‑1a hash – public‑domain implementation *
 * -------------------------------------------------- */
static uint64_t wl_hash_cache = 0;
static uint64_t fnv1a64(const void *data, size_t n_bytes)
{
    const uint8_t *p = (const uint8_t*)data;
    uint64_t hash = 0xcbf29ce484222325ULL;      /* FNV offset basis */
    while (n_bytes--) {
        hash ^= *p++;
        hash *= 0x100000001b3ULL;              /* FNV prime */
    }
    return hash;
}

// ---------- Global Memory ----------
static double* wl_cache = NULL;
static size_t cached_n_wl = 0;

static double* cached_a_w = NULL;
static double* cached_a0 = NULL;
static double* cached_a1 = NULL;
static double* cached_bb_w = NULL;
static double* cached_r_rs_b = NULL;

static double* a_w_wl = NULL, *a_w_val = NULL;
static size_t a_w_wl_n = 0;
static double* a0a1_wl = NULL, *a0_val = NULL, *a1_val = NULL;
static size_t a0a1_wl_n = 0;
static double* r_rs_b_wl = NULL, *r_rs_b_matrix = NULL;
static char** r_rs_b_class_names = NULL;
static size_t r_rs_b_class_n = 0, r_rs_b_wl_n = 0;

// ---------- Loaders ----------

int load_pure_water(const double* wl, const double* a, size_t n) {
    if (!wl || !a || n == 0) return 1;
    a_w_wl_n = n;
    a_w_wl = realloc(a_w_wl, sizeof(double) * n);
    a_w_val = realloc(a_w_val, sizeof(double) * n);
    if (!a_w_wl || !a_w_val) return 2;
    memcpy(a_w_wl, wl, sizeof(double) * n);
    memcpy(a_w_val, a, sizeof(double) * n);
    return 0;
}

int load_a0_a1(const double* wl, const double* a0, const double* a1, size_t n) {
    if (!wl || !a0 || !a1 || n == 0) return 1;
    a0a1_wl_n = n;
    a0a1_wl = realloc(a0a1_wl, sizeof(double) * n);
    a0_val = realloc(a0_val, sizeof(double) * n);
    a1_val = realloc(a1_val, sizeof(double) * n);
    if (!a0a1_wl || !a0_val || !a1_val) return 2;
    memcpy(a0a1_wl, wl, sizeof(double) * n);
    memcpy(a0_val, a0, sizeof(double) * n);
    memcpy(a1_val, a1, sizeof(double) * n);
    return 0;
}

static void free_class_names(void)
{
    if (r_rs_b_class_names) {
        for (size_t j = 0; j < r_rs_b_class_n; ++j)
            free(r_rs_b_class_names[j]);
        free(r_rs_b_class_names);
    }
    r_rs_b_class_names = NULL;
    r_rs_b_class_n     = 0;
}

int load_r_rs_b(const double  *wl,
                const char   **class_names,
                const double  *matrix,
                size_t         wl_n,
                size_t         class_n)
{
    /* 1. Allocate new blocks ------------------------------------ */
    double *tmp_wl     = realloc(r_rs_b_wl, sizeof(double) * wl_n);
    double *tmp_matrix = realloc(r_rs_b_matrix, sizeof(double) * wl_n * class_n);
    if (!tmp_wl || !tmp_matrix) return 1;

    char **tmp_names = malloc(sizeof(char*) * class_n);
    if (!tmp_names) return 1;

    for (size_t j = 0; j < class_n; ++j) {
        tmp_names[j] = strdup(class_names[j]);
        if (!tmp_names[j]) {
            for (size_t k = 0; k < j; ++k) free(tmp_names[k]);
            free(tmp_names);
            return 1;
        }
    }

    /* 2. Commit + copy data ------------------------------------- */
    free_class_names();                 /* release previous names */
    r_rs_b_wl         = tmp_wl;
    r_rs_b_matrix     = tmp_matrix;
    r_rs_b_class_names = tmp_names;
    r_rs_b_class_n = class_n;
    r_rs_b_wl_n = wl_n;

    memcpy(r_rs_b_wl,     wl, sizeof(double) * wl_n);
    memcpy(r_rs_b_matrix, matrix, sizeof(double) * wl_n * class_n);

    return 0;
}

// ---------- Interpolation Functions ----------

//double interpolate_scalar(const double* wl, const double* val, size_t n, double target) {
//    if (target < wl[0] || target > wl[n-1]) return 0.0;
//    for (size_t i = 0; i < n - 1; i++) {
//        if (target >= wl[i] && target <= wl[i+1]) {
//            double slope = (val[i+1] - val[i]) / (wl[i+1] - wl[i]);
//            return val[i] + slope * (target - wl[i]);
//        }
//    }
//    return 0.0;
//}

double interpolate_scalar(const double *wl,
                          const double *val,
                          size_t        n,
                          double        target)
{
    /* Return 0 if target is outside the table ------------------ */
    if (target < wl[0]      ) return 0.0;         /* below range  */
    if (target > wl[n-1]    ) return 0.0;         /* above range  */

    /* Exact matches at the boundaries -------------------------- */
    if (target == wl[0])     return val[0];
    if (target == wl[n-1])   return val[n-1];

    /* Binary search for enclosing interval --------------------- */
    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (target < wl[mid])
            hi = mid;
        else
            lo = mid;
    }
    /* Linear interpolation within [wl[lo], wl[hi]] */
    double t = (target - wl[lo]) / (wl[hi] - wl[lo]);
    return val[lo] + t * (val[hi] - val[lo]);

//    for (int i = 0; i < n - 1; i++) {
//        if (wl >= wl_arr[i] && wl <= wl_arr[i+1]) {
//            double slope = (val_arr[i+1] - val_arr[i]) / (wl_arr[i+1] - wl_arr[i]);
//            return val_arr[i] + slope * (wl - wl_arr[i]);
//        }
//    }

}

int interpolate_vector(const double* wl_old, const double* val_old,
                       const double* wl_new, size_t n_old, size_t n_new,
                       double* val_out) {
    if (wl_old[0] > wl_new[0] || wl_old[n_old-1] < wl_new[cached_n_wl-1])
        fprintf(stderr,
                "build_cache: target grid extends beyond source table (%g–%g vs %g–%g)\n",
                wl_new[0], wl_new[cached_n_wl-1], wl_old[0], wl_old[n_old-1]);

    for (size_t i = 0; i < n_new; i++) {
        val_out[i] = interpolate_scalar(wl_old, val_old, n_old, wl_new[i]);
    }
    return 0;
}

int interpolate_matrix(const double* wl_old, const double* mat_old,
                       const double* wl_new, size_t n_old, size_t n_new, size_t n_class,
                       double* mat_out) {
    if (wl_old[0] > wl_new[0] || wl_old[n_old-1] < wl_new[cached_n_wl-1])
        fprintf(stderr,
                "build_cache: target grid extends beyond source table (%g–%g vs %g–%g)\n",
                wl_new[0], wl_new[cached_n_wl-1], wl_old[0], wl_old[n_old-1]);

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

    double *tmp = realloc(wl_cache, sizeof(double) * n);
    if (!tmp) return 3;
    wl_cache = tmp;

    memcpy(wl_cache, wl, sizeof(double) * n);
    cached_n_wl = n;

    tmp = realloc(cached_a0, sizeof(double) * n);
    if (!tmp) return 3;
    cached_a0 = tmp;
    tmp = realloc(cached_a1, sizeof(double) * n);
    if (!tmp) return 3;
    cached_a1 = tmp;
    tmp = realloc(cached_a_w, sizeof(double) * n);
    if (!tmp) return 3;
    cached_a_w = tmp;
    tmp = realloc(cached_bb_w, sizeof(double) * n);
    if (!tmp) return 3;
    cached_bb_w = tmp;
    tmp = realloc(cached_r_rs_b, sizeof(double) * n * r_rs_b_class_n);
    if (!tmp) return 3;
    cached_r_rs_b = tmp;

    interpolate_vector(a_w_wl, a_w_val, wl, a_w_wl_n, cached_n_wl, cached_a_w);
    interpolate_vector(a0a1_wl, a0_val, wl, a0a1_wl_n, cached_n_wl, cached_a0);
    interpolate_vector(a0a1_wl, a1_val, wl, a0a1_wl_n, cached_n_wl, cached_a1);
    interpolate_matrix(r_rs_b_wl, r_rs_b_matrix, wl, r_rs_b_wl_n, cached_n_wl, r_rs_b_class_n, cached_r_rs_b);

    for (size_t i = 0; i < n; i++) {
        double lambda = wl[i];
        double b1 = 0.00111;
        double lambda1 = 500.0;
        double exponent = -4.32;
        cached_bb_w[i] = b1 * pow(lambda / lambda1, exponent);
    }

    wl_hash_cache = fnv1a64(wl, n * sizeof(double));
    return 0;
}

/*-------------------------------------------------------------*
 *  Cache guard / validation utilities                         *
 *-------------------------------------------------------------*/

/* return codes:
 *  0  – OK, cache ready (built or already compatible)
 *  1  – spectral tables not loaded yet
 *  2  – wavelength list differs from existing cache
 *  3  – build_cache() failed (propagates its error)
 */
int ensure_cache(const double *wl, size_t n)
{
    /* 1.  Are the master tables in memory? */
    if (!a0a1_wl || !a_w_wl || !r_rs_b_wl) return 1;

    if (cached_n_wl == n && wl_hash_cache) {
        uint64_t h = fnv1a64(wl, n * sizeof(double));
        if (h == wl_hash_cache)
            return 0;
    }

    /* 3.  Build (or rebuild) the cache --------------------------- */
    int rc = build_cache(wl, n);
    if (rc) return 3;

    return 0;
}

// ---------- Cache Accessors ----------

const double* get_a_w()      { return cached_a_w; }
const double* get_a0()       { return cached_a0; }
const double* get_a1()       { return cached_a1; }
const double* get_bb_w()     { return cached_bb_w; }
const double* get_r_rs_b()   { return cached_r_rs_b; }
const char**  get_r_rs_b_class_names() { return (const char**)r_rs_b_class_names; }
size_t get_n_wl()              { return cached_n_wl; }
size_t get_n_class()      { return r_rs_b_class_n; }

/* Free every dynamically allocated table so valgrind stays quiet */
void saber_reset_tables(void)
{
    /* wavelength & optics tables */
    free(a_w_wl);      free(a_w_val);
    free(a0a1_wl);     free(a0_val);  free(a1_val);
    free(r_rs_b_wl);   free(r_rs_b_matrix);

    /* class names */
    if (r_rs_b_class_names) {
        for (size_t j = 0; j < r_rs_b_class_n; ++j)
            free(r_rs_b_class_names[j]);
        free(r_rs_b_class_names);
    }

    /* cached spectra */
    free(wl_cache);   free(cached_a_w);  free(cached_a0);
    free(cached_a1);  free(cached_bb_w); free(cached_r_rs_b);

    /* reset globals */
    a_w_wl = a_w_val = a0a1_wl = a0_val = a1_val = NULL;
    r_rs_b_wl = r_rs_b_matrix = NULL;
    r_rs_b_class_names = NULL;
    wl_cache = cached_a_w = cached_a0 = cached_a1 =
    cached_bb_w = cached_r_rs_b = NULL;
    cached_n_wl = r_rs_b_class_n = 0;
    wl_hash_cache = 0;
}