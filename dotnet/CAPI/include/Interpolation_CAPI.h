/*
 * Author:
 * 2025/12/07 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "XYZ_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

LNLIB_EXPORT double LNLIB_INTERPOLATION_get_total_chord_length(
        XYZ_C* through_points,
        int through_points_count
);

LNLIB_EXPORT void    LNLIB_INTERPOLATION_get_chord_parameterization(
    XYZ_C* through_points,
    int through_points_count,
    double* out_params,
    int* out_size
);

LNLIB_EXPORT double  LNLIB_INTERPOLATION_get_centripetal_length(
    XYZ_C* through_points,
    int through_points_count);

LNLIB_EXPORT void    LNLIB_INTERPOLATION_get_centripetal_parameterization(
    XYZ_C* through_points,
    int through_points_count,
    double* out_params,
    int* out_size
);

LNLIB_EXPORT void LNLIB_INTERPOLATION_average_knot_vector(
    int degree,
    double* params,
    int params_count,
    double* out_knot_vector,
    int* out_knot_vector_size
);

LNLIB_EXPORT int LNLIB_INTERPOLATION_get_surface_mesh_parameterization(
    XYZ_C* through_points,
    int rows,
    int cols,
    double* params_u,
    int* out_params_u_size,
    double* params_v,
    int* out_params_v_size);

LNLIB_EXPORT void LNLIB_INTERPOLATION_compute_tangent_by_three_points(
    XYZ_C* through_points,
    int through_points_count,
    XYZ_C* out_tangents,
    int* out_size
);


LNLIB_EXPORT int LNLIB_INTERPOLATION_compute_tangent_by_five_points(
    XYZ_C* through_points,
    int through_points_count,
    XYZ_C* out_tangents,
    int* out_size
);


LNLIB_EXPORT int LNLIB_INTERPOLATION_computer_weight_for_rational_quadratic_interpolation(
    XYZ_C start_point,
    XYZ_C middle_control_point,
    XYZ_C end_point,
    double* out_weight
);

LNLIB_EXPORT void LNLIB_INTERPOLATION_compute_knot_vector(
    int degree,
    int control_points_count,
    double* params,
    int params_count,
    double* out_knot_vector,
    int* out_knot_vector_size
);


#ifdef __cplusplus
}
#endif