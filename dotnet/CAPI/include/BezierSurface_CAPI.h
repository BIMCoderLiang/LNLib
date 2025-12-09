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
#include "UV_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

LNLIB_EXPORT XYZ_C LNLIB_BEZIERSURF_get_point_on_surface_by_deCasteljau(
    int degree_u,
    int degree_v,
    const XYZ_C* control_points,
    int cp_rows,
    int cp_cols,
    UV_C uv);

LNLIB_EXPORT XYZW_C LNLIB_BEZIERSURF_get_rational_point_on_surface_by_deCasteljau(
    int degree_u,
    int degree_v,
    const XYZW_C* control_points,
    int cp_rows,
    int cp_cols,
    UV_C uv);

#ifdef __cplusplus
}
#endif