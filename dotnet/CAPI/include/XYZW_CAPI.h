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
#include "XYZ_CAPI.h"
#include "XYZW.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct { double wx, wy, wz, w; } XYZW_C;

LNLIB_EXPORT XYZW_C LNLIB_XYZW_create(double wx, double wy, double wz, double w);
LNLIB_EXPORT XYZW_C LNLIB_XYZW_create_from_xyz(XYZ_C xyz, double w);

LNLIB_EXPORT double LNLIB_XYZW_get_wx(XYZW_C xyzw);
LNLIB_EXPORT double LNLIB_XYZW_get_wy(XYZW_C xyzw);
LNLIB_EXPORT double LNLIB_XYZW_get_wz(XYZW_C xyzw);
LNLIB_EXPORT double LNLIB_XYZW_get_w(XYZW_C xyzw);
LNLIB_EXPORT XYZW_C LNLIB_XYZW_set_w(XYZW_C xyzw, double w);

LNLIB_EXPORT XYZ_C	LNLIB_XYZW_to_xyz(XYZW_C xyzw, int divideWeight);
LNLIB_EXPORT int	LNLIB_XYZW_is_almost_equal(XYZW_C a, XYZW_C b, double epsilon);
LNLIB_EXPORT double LNLIB_XYZW_distance(XYZW_C a, XYZW_C b);

LNLIB_EXPORT XYZW_C LNLIB_XYZW_add(XYZW_C a, XYZW_C b);
LNLIB_EXPORT XYZW_C LNLIB_XYZW_subtract(XYZW_C a, XYZW_C b);
LNLIB_EXPORT XYZW_C LNLIB_XYZW_multiply(XYZW_C a, double scalar);
LNLIB_EXPORT XYZW_C LNLIB_XYZW_divide(XYZW_C a, double scalar);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
inline LNLib::XYZW ToXYZW(XYZW_C c) { return LNLib::XYZW(c.wx, c.wy, c.wz, c.w); }
inline XYZW_C FromXYZW(const LNLib::XYZW& v) { return { v.GetWX(), v.GetWY(), v.GetWZ(), v.GetW() }; }
#endif