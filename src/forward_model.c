#include "forward_model.h"
#include "snell_law.h"
#include <math.h>

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
        } else if (water_type == 2) {
            f_rs = 0.0512 *
                   (1 + 4.6659 * omega_b +
                    -7.8387 * omega_b * omega_b +
                    5.4571 * omega_b * omega_b * omega_b) *
                   (1 + (0.1098 / cos(sun_w_rad))) *
                   (1 + (0.4021 / cos(view_w_rad)));
        } else {
            return 3;
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

/*-------------------------------------------------------------------------*/
/*  Recover bottom reflectance r_b(λ) from observed Rrs(λ)                 */
/*                                                                         */
/*  Returns 0 on success, >0 on error:                                     */
/*      1 – null pointer                                                   */
/*      2 – invalid geometry/parameters (h_w ≤ 0 or deep water requested)  */
/*      3 – denominator ≈ 0 (numerically unstable)                         */
/*-------------------------------------------------------------------------*/
int retrieve_r_rs_b_am03(
        const double *wavelength,   /* [n] λ  (nm)               */
        const double *a,            /* [n] absorption a(λ)       */
        const double *bb,           /* [n] backscatter bb(λ)     */
        const double *r_rs_obs,      /* [n] observed Rrs(λ)       */
        size_t        n,            /* number of bands           */
        int           water_type,   /* 0 = clear, 1 = turbid     */
        double        theta_sun_deg,
        double        theta_view_deg,
        double        h_w,          /* water depth  (m)          */
        double       *r_rs_b_out       /* [n]  ← recovered r_b(λ)   */
)
{
    if (!wavelength || !a || !bb || !r_rs_obs || !r_rs_b_out) return 1;
    if (h_w <= 0.0) return 2;           /* bottom retrieval only makes sense for shallow water */

    /* --- 1. Snell conversion of angles to water column ----------------- */
    double view_w_rad = 0.0, sun_w_rad = 0.0;
    snell_law(theta_view_deg, theta_sun_deg, &view_w_rad, &sun_w_rad);

    /* --- 2. Constant coefficients (Albert & Mobley 2003) --------------- */
    const double Ars1 = 1.1576;
    const double Ars2 = 1.0389;
    const double k0   = (water_type == 1) ? 1.0395 : 1.0546;

    for (size_t i = 0; i < n; ++i) {

        const double ext = a[i] + bb[i];        /* total attenuation a+bb */
        if (ext <= 0.0) {                       /* avoid division by zero  */
            r_rs_b_out[i] = 0.0;
            continue;
        }

        const double omega_b = bb[i] / ext;     /* single-backscatter albedo */

        /* ----- 2a. Fresnel/geometry factor f_rs (deep water) ----------- */
        double f_rs;
        if (water_type == 1) {                  /* turbid case */
            f_rs = 0.095;
        } else if (water_type == 2) {                                /* clear / case-1 water  */
            f_rs = 0.0512 *
                   (1 + 4.6659  * omega_b
                    - 7.8387 * omega_b * omega_b
                    + 5.4571 * pow(omega_b, 3.0)) *
                   (1 + 0.1098 / cos(sun_w_rad)) *
                   (1 + 0.4021 / cos(view_w_rad));
        } else {
            return 3;
        }

        const double rrs_deep = f_rs * omega_b;

        /* ----- 2b. Diffuse attenuation coefficients ------------------- */
        const double Kd  = k0 * (ext / cos(sun_w_rad));

        const double kuW = (ext / cos(view_w_rad))
                           * pow(1 + omega_b, 3.5421)
                           * (1 - 0.2786 / cos(sun_w_rad));

        const double kuB = (ext / cos(view_w_rad))
                           * pow(1 + omega_b, 2.2658)
                           * (1 - 0.0577 / cos(sun_w_rad));

        /* Pre-compute the exponential terms once */
        const double exp_W = exp(-h_w * (Kd + kuW));
        const double exp_B = exp(-h_w * (Kd + kuB));

        /* ----- 3. Invert the shallow-water equation for r_b ------------ */
        const double numerator   = r_rs_obs[i] - rrs_deep * (1.0 - Ars1 * exp_W);
        const double denominator = Ars2 * exp_B;

        if (fabs(denominator) < 1e-12) {        /* guard against blow-ups */
            r_rs_b_out[i] = 0.0;
            return 4;
        }

        r_rs_b_out[i] = numerator / denominator;

    }

    return 0;
}
