#include "../include/data_cache.h"
#include "../include/snell_law.h"
#include "../include/forward_model.h"
#include "../include/iop_from_oac.h"
#include "../include/r_rs_b_lmm.h"

#include <stdio.h>
#include <string.h>

// A simple helper to print a vector
void print_vector(const char* label, const double* vec, size_t n) {
    printf("%s = [", label);
    for (size_t i = 0; i < n; i++) {
        printf("%.4f", vec[i]);
        if (i < n - 1) printf(", ");
    }
    printf("]\n");
}

int main() {
    size_t n = 3;
    double wl[] = {443.0, 490.0, 555.0};

    // 1. Set up synthetic absorption data
    double a_w[] = {0.01, 0.02, 0.03};
    double a0[] = {0.06, 0.05, 0.04};
    double a1[] = {0.01, 0.01, 0.01};
    load_pure_water(wl, a_w, n);
    load_a0_a1(wl, a0, a1, n);

    // 2. Set up synthetic bottom reflectance matrix
    const char *class_names[] = {"sand", "algae"};
    size_t n_cls = 2;
    double rrs_b[] = {
            0.01, 0.03,
            0.02, 0.04,
            0.03, 0.05
    };
    load_bottom_reflectance(wl, rrs_b, class_names, n, n_cls);

    // 3. Build cache
    build_cache(wl, n);

    // 4. Run IOP model
    const char *pnames[] = {"chl", "a_g_440", "bb_p_550"};
    double pvals[] = {1.0, 0.1, 0.01};
    double a_out[n], bb_out[n];
    iop_from_oac(wl, n, pnames, pvals, 3, a_out, bb_out);
    print_vector("a", a_out, n);
    print_vector("bb", bb_out, n);

    // 5. Compute r_rs_b from fractions
    const char *mix_classes[] = {"sand", "algae"};
    double mix_fractions[] = {0.7, 0.3};
    double r_rs_b[n];
    compute_r_rs_b_lmm(mix_classes, mix_fractions, 2, r_rs_b);
    print_vector("r_rs_b", r_rs_b, n);

    // 6. Run forward model
    double rrs_0m[3];
    forward_am03(wl, a_out, bb_out, n,
                       2, 20.0, 0.0, 1, 5.0, r_rs_b, rrs_0m);
    print_vector("Rrs", rrs_0m, n);

    return 0;
}