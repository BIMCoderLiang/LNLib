/*
 * Author:
 * 2025/11/25 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once

#include "LNLibDefinitions.h"
#include "UV.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double u;
    double v;
} UV_C;

LNLIB_EXPORT UV_C   LNLIB_UV_create(double u, double v);
LNLIB_EXPORT double LNLIB_UV_get_u(UV_C uv);
LNLIB_EXPORT double LNLIB_UV_get_v(UV_C uv);
LNLIB_EXPORT UV_C   LNLIB_UV_set_u(UV_C uv, double u);
LNLIB_EXPORT UV_C   LNLIB_UV_set_v(UV_C uv, double v);

LNLIB_EXPORT int    LNLIB_UV_is_zero(UV_C uv, double epsilon);
LNLIB_EXPORT int    LNLIB_UV_is_unit(UV_C uv, double epsilon);
LNLIB_EXPORT int    LNLIB_UV_is_almost_equal(UV_C a, UV_C b, double epsilon);
LNLIB_EXPORT double LNLIB_UV_length(UV_C uv);
LNLIB_EXPORT double LNLIB_UV_sqr_length(UV_C uv);
LNLIB_EXPORT UV_C   LNLIB_UV_normalize(UV_C uv);
LNLIB_EXPORT UV_C   LNLIB_UV_add(UV_C a, UV_C b);
LNLIB_EXPORT UV_C   LNLIB_UV_subtract(UV_C a, UV_C b);
LNLIB_EXPORT UV_C   LNLIB_UV_negative(UV_C uv);
LNLIB_EXPORT double LNLIB_UV_dot_product(UV_C a, UV_C b);
LNLIB_EXPORT double LNLIB_UV_cross_product(UV_C a, UV_C b);
LNLIB_EXPORT double LNLIB_UV_distance(UV_C a, UV_C b);

LNLIB_EXPORT UV_C   LNLIB_UV_scale(UV_C uv, double factor);
LNLIB_EXPORT UV_C   LNLIB_UV_divide(UV_C uv, double divisor);


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
inline LNLib::UV ToUV(UV_C v) { return LNLib::UV(v.u, v.v); }
inline UV_C FromUV(const LNLib::UV& v) { return { v.U(), v.V() }; }
#endif