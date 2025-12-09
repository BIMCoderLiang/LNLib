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
#include "LNEnums_CAPI.h"
#include "XYZ_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

LNLIB_EXPORT LNLIB_ENUMS_CurveCurveIntersectionType_C LNLIB_INTERSECT_compute_rays(
    XYZ_C point_0, XYZ_C vector_0,
    XYZ_C point_1, XYZ_C vector_1,
    double* param_0,
    double* param_1,
    XYZ_C* intersect_point);


LNLIB_EXPORT LNLIB_ENUMS_LinePlaneIntersectionType_C LNLIB_INTERSECT_compute_line_and_plane(
    XYZ_C normal,
    XYZ_C point_on_plane,
    XYZ_C point_on_line,
    XYZ_C line_direction,
    XYZ_C* intersect_point);

#ifdef __cplusplus
}
#endif