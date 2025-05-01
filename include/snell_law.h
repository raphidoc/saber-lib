//
// Created by raphael on 5/1/25.
//

#ifndef SABER_LIB_SNELL_LAW_H
#define SABER_LIB_SNELL_LAW_H

#ifdef __cplusplus
extern "C" {
#endif

void snell_law(double theta_view_deg, double theta_sun_deg,
                     double* view_w, double* sun_w);

#ifdef __cplusplus
}
#endif

#endif //SABER_LIB_SNELL_LAW_H
