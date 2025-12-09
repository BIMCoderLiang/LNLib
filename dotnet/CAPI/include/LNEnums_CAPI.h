/*
 * Author:
 * 2025/11/29 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

 #pragma once

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    CURVE_CURVE_INTERSECTING = 0,
    CURVE_CURVE_PARALLEL = 1,
    CURVE_CURVE_COINCIDENT = 2,
    CURVE_CURVE_SKEW = 3
} LNLIB_ENUMS_CurveCurveIntersectionType_C;

typedef enum {
    LINE_PLANE_INTERSECTING = 0,
    LINE_PLANE_PARALLEL = 1,
    LINE_PLANE_ON = 2
} LNLIB_ENUMS_LinePlaneIntersectionType_C;

typedef enum {
    CURVE_NORMAL_NORMAL = 0,
    CURVE_NORMAL_BINORMAL = 1
} LNLIB_ENUMS_CurveNormal_C;

typedef enum {
    SURFACE_DIRECTION_ALL = 0,
    SURFACE_DIRECTION_U = 1,
    SURFACE_DIRECTION_V = 2
} LNLIB_ENUMS_SurfaceDirection_C;

typedef enum {
    SURFACE_CURVATURE_MAXIMUM = 0,
    SURFACE_CURVATURE_MINIMUM = 1,
    SURFACE_CURVATURE_GAUSS = 2,
    SURFACE_CURVATURE_MEAN = 3,
    SURFACE_CURVATURE_ABS = 4,
    SURFACE_CURVATURE_RMS = 5
} LNLIB_ENUMS_SurfaceCurvature_C;

typedef enum {
    INTEGRATOR_SIMPSON = 0,
    INTEGRATOR_GAUSS_LEGENDRE = 1,
    INTEGRATOR_CHEBYSHEV = 2
} LNLIB_ENUMS_IntegratorType_C;

typedef enum {
    OFFSET_TILLER_HANSON = 0,
    OFFSET_PIEGL_TILLER = 1
} LNLIB_ENUMS_OffsetType_C;

#ifdef __cplusplus
}
#endif