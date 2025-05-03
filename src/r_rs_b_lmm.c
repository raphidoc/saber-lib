#include "r_rs_b_lmm.h"
#include "data_cache.h"
#include <string.h>
#include <stdio.h>

int compute_r_rs_b_lmm(
        const char** class_names, const double* class_fractions, size_t n_frac,
        double* out_r_rs_b
) {
    if (!class_names || !class_fractions || !out_r_rs_b) return 1;

    const double* r_rs_b = get_r_rs_b();
    const char** colnames = get_r_rs_b_class_names();
    size_t n_wl = get_n_wl();
    size_t n_class = get_n_class();

    if (!r_rs_b || !colnames || n_wl == 0 || n_class == 0) return 2;

    // Zero initialize
    for (size_t i = 0; i < n_wl; i++) out_r_rs_b[i] = 0.0;

    // Normalize fractions
    double sum = 0.0;
    for (size_t j = 0; j < n_frac; j++) {
        if (class_fractions[j] < 0.0) return 3;  // invalid
        sum += class_fractions[j];
    }

    if (sum <= 0.0) return 4;

    for (size_t j = 0; j < n_frac; j++) {
        const char* name = class_names[j];
        double weight = class_fractions[j] / sum;

        int matched = -1;
        for (size_t k = 0; k < n_class; k++) {
            if (strcmp(name, colnames[k]) == 0) {
                matched = (int)k;
                break;
            }
        }

        if (matched < 0) {
            fprintf(stderr, "Class name '%s' not found in cached bottom reflectance\n", name);
            return 5;
        }

        // Add weighted column to output
        for (size_t i = 0; i < n_wl; i++) {
            out_r_rs_b[i] += weight * r_rs_b[i + matched * n_wl];
        }
    }

    return 0;
}
