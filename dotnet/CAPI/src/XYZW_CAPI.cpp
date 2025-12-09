/*
 * Author:
 * 2025/11/26 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "XYZW_CAPI.h"
#include "XYZW.h"
#include "XYZ.h"
#include "MathUtils.h"


XYZW_C LNLIB_XYZW_create(double wx, double wy, double wz, double w)
{
	return { wx, wy, wz, w };
}

XYZW_C LNLIB_XYZW_create_from_xyz(XYZ_C xyz, double w)
{
	return { w * xyz.x, w * xyz.y, w * xyz.z, w };
}

XYZ_C LNLIB_XYZW_to_xyz(XYZW_C v, int divideWeight)
{
	LNLib::XYZW vw(v.wx, v.wy, v.wz, v.w);
	LNLib::XYZ res = vw.ToXYZ(divideWeight);
	return { res.GetX(), res.GetY(), res.GetZ() };
}

int LNLIB_XYZW_is_almost_equal(XYZW_C a, XYZW_C b, double epsilon)
{
	return
		(LNLib::MathUtils::IsAlmostEqualTo(a.wx, b.wx, epsilon) &&
			LNLib::MathUtils::IsAlmostEqualTo(a.wy, b.wy, epsilon) &&
			LNLib::MathUtils::IsAlmostEqualTo(a.wz, b.wz, epsilon) &&
			LNLib::MathUtils::IsAlmostEqualTo(a.w, b.w, epsilon)) ? 1 : 0;
}

XYZW_C LNLIB_XYZW_add(XYZW_C a, XYZW_C b)
{
	LNLib::XYZW va(a.wx, a.wy, a.wz, a.w);
	LNLib::XYZW vb(b.wx, b.wy, b.wz, b.w);
	LNLib::XYZW res = va + vb;
	return { res.GetWX(), res.GetWY(), res.GetWZ(), res.GetW() };
}

XYZW_C LNLIB_XYZW_subtract(XYZW_C a, XYZW_C b)
{
	LNLib::XYZW va(a.wx, a.wy, a.wz, a.w);
	LNLib::XYZW vb(b.wx, b.wy, b.wz, b.w);
	LNLib::XYZW res = va - vb;
	return { res.GetWX(), res.GetWY(), res.GetWZ(), res.GetW() };
}

XYZW_C LNLIB_XYZW_multiply(XYZW_C xyzw, double scalar)
{
	LNLib::XYZW va(xyzw.wx, xyzw.wy, xyzw.wz, xyzw.w);
	LNLib::XYZW res = va * scalar;
	return { res.GetWX(), res.GetWY(), res.GetWZ(), res.GetW() };
}

XYZW_C LNLIB_XYZW_divide(XYZW_C xyzw, double scalar)
{
	LNLib::XYZW va(xyzw.wx, xyzw.wy, xyzw.wz, xyzw.w);
	LNLib::XYZW res = va / scalar;
	return { res.GetWX(), res.GetWY(), res.GetWZ(), res.GetW() };
}

double LNLIB_XYZW_distance(XYZW_C a, XYZW_C b)
{
	LNLib::XYZW va(a.wx, a.wy, a.wz, a.w);
	LNLib::XYZW vb(b.wx, b.wy, b.wz, b.w);
	return va.Distance(vb);
}

double LNLIB_XYZW_get_wx(XYZW_C xyzw) { return xyzw.wx; }
double LNLIB_XYZW_get_wy(XYZW_C xyzw) { return xyzw.wy; }
double LNLIB_XYZW_get_wz(XYZW_C xyzw) { return xyzw.wz; }
double LNLIB_XYZW_get_w(XYZW_C xyzw) { return xyzw.w; }

XYZW_C LNLIB_XYZW_set_w(XYZW_C xyzw, double w)
{
	LNLib::XYZW vw(xyzw.wx, xyzw.wy, xyzw.wz, xyzw.w);
	LNLib::XYZ res = vw.ToXYZ(1);
	return { w * res.GetX(), w * res.GetY(), w * res.GetZ(), w };
}