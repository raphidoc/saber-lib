#include "../include/snell_law.h"
#include <math.h>

// Cached last-used input/output
static double cached_theta_view = -9999;
static double cached_theta_sun  = -9999;
static double cached_view_w     = 0;
static double cached_sun_w      = 0;

void snell_law(double theta_view_deg, double theta_sun_deg,
                     double* view_w, double* sun_w) {
    if (theta_view_deg == cached_theta_view && theta_sun_deg == cached_theta_sun) {
        *view_w = cached_view_w;
        *sun_w  = cached_sun_w;
        return;
    }

    double theta_view_rad = theta_view_deg * M_PI / 180.0;
    double theta_sun_rad  = theta_sun_deg  * M_PI / 180.0;

    double n_air = 1.0;
    double n_water = 1.33;

    *view_w = asin((n_air / n_water) * sin(theta_view_rad));
    *sun_w  = asin((n_air / n_water) * sin(theta_sun_rad));

    // Cache result
    cached_theta_view = theta_view_deg;
    cached_theta_sun  = theta_sun_deg;
    cached_view_w     = *view_w;
    cached_sun_w      = *sun_w;
}
