/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "XYZW_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

LNLIB_EXPORT int LNLIB_VALID_is_valid_knotVector(const double* knot_vector, int knot_vector_count);
LNLIB_EXPORT int LNLIB_VALID_is_valid_bezier(int degree, int control_points_count);
LNLIB_EXPORT int LNLIB_VALID_is_valid_bspline(int degree, int knot_vector_count, int control_points_count);
LNLIB_EXPORT int LNLIB_VALID_is_valid_nurbs(int degree, int knot_vector_count, int control_points_count);

LNLIB_EXPORT double LNLIB_VALID_compute_curve_modify_tolerance(
    const XYZW_C* control_points,
    int control_points_count
);

#ifdef __cplusplus
}
#endif