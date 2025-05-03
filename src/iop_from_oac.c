#include "iop_from_oac.h"
#include "data_cache.h"
#include <string.h>
#include <math.h>

static double get_named_value(const char* name, const char** names, const double* values, size_t n, int* found) {
    for (size_t i = 0; i < n; i++) {
        if (strcmp(names[i], name) == 0) {
            *found = 1;
            return values[i];
        }
    }
    *found = 0;
    return 0.0;
}

int iop_from_oac(
        const double* wavelength, size_t n,
        const char** param_names, const double* param_values, size_t n_param,
        double* a_out, double* bb_out
) {
    if (!wavelength || !a_out || !bb_out) return 1;
    int rc = ensure_cache(wavelength, n);
    if (rc) {
        return rc;
    }

    const double* aw_ptr   = get_a_w();
    const double* a0_ptr   = get_a0();
    const double* a1_ptr   = get_a1();
    const double* bb_w_ptr = get_bb_w();

    // Fetch named parameters
    int found;
    double chl          = get_named_value("chl", param_names, param_values, n_param, &found);
    int has_chl         = found;

    double a_g_440      = get_named_value("a_g_440", param_names, param_values, n_param, &found);
    int has_a_g_440     = found;

    double a_nap_440    = get_named_value("a_nap_440", param_names, param_values, n_param, &found);
    int has_a_nap_440   = found;

    double bb_p_550     = get_named_value("bb_p_550", param_names, param_values, n_param, &found);
    int has_bb_p_550    = found;

    double a_g_s_g      = get_named_value("a_g_s_g", param_names, param_values, n_param, &found);
    double a_g_s_d      = get_named_value("a_g_s_d", param_names, param_values, n_param, &found);
    int has_a_g_slopes  = found;

    double a_nap_s_d    = get_named_value("a_nap_s_d", param_names, param_values, n_param, &found);
    int has_a_nap_slope = found;

    double bb_p_gamma   = get_named_value("bb_p_gamma", param_names, param_values, n_param, &found);
    int has_bb_p_gamma  = found;

    // Compute loop
    for (size_t i = 0; i < n; i++) {
        double wl = wavelength[i];

        // Phytoplankton absorption
        double a_phy = 0.0;
        if (has_chl) {
            double aph_440 = 0.06 * pow(chl, 0.65);
            double a0 = a0_ptr[i];
            double a1 = a1_ptr[i];
            a_phy = (a0 + a1 * log(aph_440)) * aph_440;
            if (a_phy < 0.0) a_phy = 0.0;
        }

        // CDOM absorption
        double a_g = 0.0;
        if (has_a_g_440) {
            double slope = has_a_g_slopes ? (a_g_s_g + a_g_s_d) : 0.017;
            a_g = a_g_440 * exp(-slope * (wl - 440.0));
        }

        // NAP absorption
        double a_nap = 0.0;
        if (has_a_nap_440) {
            double slope = has_a_nap_slope ? a_nap_s_d : 0.0116;
            a_nap = a_nap_440 * exp(-slope * (wl - 440.0));
        }

        // Particle backscattering
        double bb_p = 0.0;
        if (has_bb_p_550) {
            double gamma = has_bb_p_gamma ? bb_p_gamma : 0.46;
            bb_p = bb_p_550 * pow(wl / 550.0, -gamma);
        }

        a_out[i]  = aw_ptr[i] + a_phy + a_g + a_nap;
        bb_out[i] = bb_w_ptr[i] + bb_p;
    }

    return 0;
}
