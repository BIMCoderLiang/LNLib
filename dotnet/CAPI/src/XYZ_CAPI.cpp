/*
 * Author:
 * 2025/11/25 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "XYZ_CAPI.h"
#include "XYZ.h"
#include "MathUtils.h"
#include <cmath>

extern "C" {

	LNLIB_EXPORT XYZ_C xyz_create(double x, double y, double z)
	{
		return { x, y, z };
	}

	LNLIB_EXPORT XYZ_C xyz_zero()
	{
		return { 0.0, 0.0, 0.0 };
	}

	LNLIB_EXPORT XYZ_C xyz_add(XYZ_C a, XYZ_C b)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ vb(b.x, b.y, b.z);
		LNLib::XYZ res = va + vb;
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT XYZ_C xyz_subtract(XYZ_C a, XYZ_C b)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ vb(b.x, b.y, b.z);
		LNLib::XYZ res = va - vb;
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT XYZ_C xyz_negative(XYZ_C a)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ res = -va;
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT XYZ_C xyz_multiply(XYZ_C a, double scalar)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ res = va * scalar;
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT XYZ_C xyz_divide(XYZ_C a, double scalar)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ res = va / scalar;
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT double xyz_length(XYZ_C v)
	{
		LNLib::XYZ vv(v.x, v.y, v.z);
		return vv.Length();
	}

	LNLIB_EXPORT double xyz_sqr_length(XYZ_C v)
	{
		LNLib::XYZ vv(v.x, v.y, v.z);
		return vv.SqrLength();
	}

	LNLIB_EXPORT int xyz_is_zero(XYZ_C v, double epsilon)
	{
		LNLib::XYZ vv(v.x, v.y, v.z);
		return vv.IsZero(epsilon) ? 1 : 0;
	}

	LNLIB_EXPORT int xyz_is_unit(XYZ_C v, double epsilon)
	{
		LNLib::XYZ vv(v.x, v.y, v.z);
		return vv.IsUnit(epsilon) ? 1 : 0;
	}

	LNLIB_EXPORT XYZ_C xyz_normalize(XYZ_C v)
	{
		LNLib::XYZ vv(v.x, v.y, v.z);
		LNLib::XYZ res = vv.Normalize();
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT double xyz_dot(XYZ_C a, XYZ_C b)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ vb(b.x, b.y, b.z);
		return va.DotProduct(vb);
	}

	LNLIB_EXPORT XYZ_C xyz_cross(XYZ_C a, XYZ_C b)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ vb(b.x, b.y, b.z);
		LNLib::XYZ res = va.CrossProduct(vb);
		return { res.GetX(), res.GetY(), res.GetZ() };
	}

	LNLIB_EXPORT double xyz_distance(XYZ_C a, XYZ_C b)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ vb(b.x, b.y, b.z);
		return va.Distance(vb);
	}

	LNLIB_EXPORT int xyz_equals(XYZ_C a, XYZ_C b)
	{
		LNLib::XYZ va(a.x, a.y, a.z);
		LNLib::XYZ vb(b.x, b.y, b.z);
		return va.IsAlmostEqualTo(vb) ? 1 : 0;
	}

}