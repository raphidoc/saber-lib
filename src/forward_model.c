#include "../include/forward_model.h"
#include "../include/snell_law.h"
#include <math.h>

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
) {
    if (!wavelength || !a || !bb || !rrs_out) return 1;
    if (shallow && (!r_b || h_w < 0)) return 2;

    // Compute viewing geometry
    double view_w_rad = 0, sun_w_rad = 0;
    snell_law(theta_view_deg, theta_sun_deg, &view_w_rad, &sun_w_rad);

    for (size_t i = 0; i < n; i++) {
        double ext = a[i] + bb[i];
        if (ext == 0) {
            rrs_out[i] = 0;
            continue;
        }

        double omega_b = bb[i] / ext;

        // Fresnel & geometry factor
        double f_rs;
        if (water_type == 1) {
            f_rs = 0.095;
        } else {
            f_rs = 0.0512 *
                   (1 + 4.6659 * omega_b +
                    -7.8387 * omega_b * omega_b +
                    5.4571 * omega_b * omega_b * omega_b) *
                   (1 + (0.1098 / cos(sun_w_rad))) *
                   (1 + (0.4021 / cos(view_w_rad)));
        }

        double rrs_deep = f_rs * omega_b;

        if (shallow) {
            double k0 = (water_type == 1) ? 1.0395 : 1.0546;

            double Kd = k0 * (ext / cos(sun_w_rad));
            double kuW = (ext / cos(view_w_rad)) *
                         pow(1 + omega_b, 3.5421) *
                         (1 - 0.2786 / cos(sun_w_rad));
            double kuB = (ext / cos(view_w_rad)) *
                         pow(1 + omega_b, 2.2658) *
                         (1 - 0.0577 / cos(sun_w_rad));

            double Ars1 = 1.1576;
            double Ars2 = 1.0389;

            rrs_out[i] = rrs_deep * (1 - (Ars1 * exp(-h_w * (Kd + kuW)))) +
                         Ars2 * r_b[i] * exp(-h_w * (Kd + kuB));
        } else {
            rrs_out[i] = rrs_deep;
        }
    }

    return 0;
}
