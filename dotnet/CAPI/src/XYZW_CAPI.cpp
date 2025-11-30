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

extern "C" {

	LNLIB_EXPORT XYZW_C xyzw_create(double wx, double wy, double wz, double w)
	{
		return { wx, wy, wz, w };
	}

	LNLIB_EXPORT XYZW_C xyzw_create_from_xyz(XYZ_C xyz, double w)
	{
		return { w * xyz.x, w * xyz.y, w * xyz.z, w };
	}

	LNLIB_EXPORT XYZ_C xyzw_to_xyz(XYZW_C v, int divideWeight)
	{
		LNLib::XYZW vw(v.wx, v.wy, v.wz, v.w);
		LNLib::XYZ res = vw.ToXYZ(divideWeight);
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT XYZW_C xyzw_add(XYZW_C a, XYZW_C b)
	{
		LNLib::XYZW va(a.wx, a.wy, a.wz, a.w);
		LNLib::XYZW vb(b.wx, b.wy, b.wz, b.w);
		LNLib::XYZW res = va + vb;
		return { res.GetWX(), res.GetWY(), res.GetWZ(), res.GetW() };
	}

	LNLIB_EXPORT XYZW_C xyzw_multiply(XYZW_C a, double scalar)
	{
		LNLib::XYZW va(a.wx, a.wy, a.wz, a.w);
		LNLib::XYZW res = va * scalar;
		return { res.GetWX(), res.GetWY(), res.GetWZ(), res.GetW() };
	}

	LNLIB_EXPORT XYZW_C xyzw_divide(XYZW_C a, double scalar)
	{
		LNLib::XYZW va(a.wx, a.wy, a.wz, a.w);
		LNLib::XYZW res = va / scalar;
		return { res.GetWX(), res.GetWY(), res.GetWZ(), res.GetW() };
	}

	LNLIB_EXPORT double xyzw_distance(XYZW_C a, XYZW_C b)
	{
		LNLib::XYZW va(a.wx, a.wy, a.wz, a.w);
		LNLib::XYZW vb(b.wx, b.wy, b.wz, b.w);
		return va.Distance(vb);
	}

	LNLIB_EXPORT double xyzw_get_wx(XYZW_C v) { return v.wx; }
	LNLIB_EXPORT double xyzw_get_wy(XYZW_C v) { return v.wy; }
	LNLIB_EXPORT double xyzw_get_wz(XYZW_C v) { return v.wz; }
	LNLIB_EXPORT double xyzw_get_w(XYZW_C v) { return v.w; }

}