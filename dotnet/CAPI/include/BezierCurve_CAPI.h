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
#include "XYZ_CAPI.h"
#include "XYZW_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

LNLIB_EXPORT XYZ_C LNLIB_BEZIERCUR_get_point_on_curve_by_bernstein(
    int degree,
    const XYZ_C* control_points,
    int control_points_count,
    double param_t);

LNLIB_EXPORT XYZW_C LNLIB_BEZIERCUR_get_rational_point_on_curve_by_bernstein(
    int degree,
    const XYZW_C* control_points,
    int control_points_count,
    double param_t);

LNLIB_EXPORT XYZ_C LNLIB_BEZIERCUR_get_point_on_curve_by_deCasteljau(
    int degree,
    const XYZ_C* control_points,
    int control_points_count,
    double param_t);

LNLIB_EXPORT XYZW_C LNLIB_BEZIERCUR_get_rational_point_on_curve_by_deCasteljau(
    int degree,
    const XYZW_C* control_points,
    int control_points_count,
    double param_t);

#ifdef __cplusplus
}
#endif