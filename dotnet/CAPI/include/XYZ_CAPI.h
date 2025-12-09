/*
 * Author:
 * 2025/11/26 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "XYZ.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double x, y, z; } XYZ_C;

LNLIB_EXPORT XYZ_C	LNLIB_XYZ_create(double x, double y, double z);

LNLIB_EXPORT double LNLIB_XYZ_get_x(XYZ_C xyz);
LNLIB_EXPORT double LNLIB_XYZ_get_y(XYZ_C xyz);
LNLIB_EXPORT double LNLIB_XYZ_get_z(XYZ_C xyz);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_set_x(XYZ_C xyz, double x);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_set_y(XYZ_C xyz, double y);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_set_z(XYZ_C xyz, double z);

LNLIB_EXPORT int	LNLIB_XYZ_is_zero(XYZ_C xyz, double epsilon);
LNLIB_EXPORT int	LNLIB_XYZ_is_unit(XYZ_C xyz, double epsilon);
LNLIB_EXPORT int	LNLIB_XYZ_is_almost_equal(XYZ_C a, XYZ_C b, double epsilon);
LNLIB_EXPORT double LNLIB_XYZ_length(XYZ_C xyz);
LNLIB_EXPORT double LNLIB_XYZ_sqr_length(XYZ_C xyz);
LNLIB_EXPORT double LNLIB_XYZ_angle_to(XYZ_C a, XYZ_C b, double epsilon);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_normalize(XYZ_C xyz);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_add(XYZ_C a, XYZ_C b);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_subtract(XYZ_C a, XYZ_C b);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_negative(XYZ_C xyz);
LNLIB_EXPORT double LNLIB_XYZ_dot_product(XYZ_C a, XYZ_C b);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_cross_product(XYZ_C a, XYZ_C b);
LNLIB_EXPORT double LNLIB_XYZ_distance(XYZ_C a, XYZ_C b);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_multiply(XYZ_C xyz, double scalar);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_divide(XYZ_C xyz, double scalar);
LNLIB_EXPORT XYZ_C	LNLIB_XYZ_create_random_orthogonal(XYZ_C xyz);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

inline LNLib::XYZ ToXYZ(XYZ_C c) { return LNLib::XYZ(c.x, c.y, c.z); }
inline XYZ_C FromXYZ(const LNLib::XYZ& v) { return { v.GetX(), v.GetY(), v.GetZ() }; }

#endif
