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


XYZ_C LNLIB_XYZ_create(double x, double y, double z)
{
	return { x, y, z };
}

double LNLIB_XYZ_get_x(XYZ_C xyz)
{
	return xyz.x;
}

double LNLIB_XYZ_get_y(XYZ_C xyz)
{
	return xyz.y;
}

double LNLIB_XYZ_get_z(XYZ_C xyz)
{
	return xyz.z;
}

XYZ_C LNLIB_XYZ_set_x(XYZ_C xyz, double x)
{
	return { x, xyz.y, xyz.z };
}

XYZ_C LNLIB_XYZ_set_y(XYZ_C xyz, double y)
{
	return { xyz.x, y, xyz.z };
}

XYZ_C LNLIB_XYZ_set_z(XYZ_C xyz, double z)
{
	return { xyz.x, xyz.y, z };
}

XYZ_C LNLIB_XYZ_add(XYZ_C a, XYZ_C b)
{
	LNLib::XYZ va(a.x, a.y, a.z);
	LNLib::XYZ vb(b.x, b.y, b.z);
	LNLib::XYZ res = va + vb;
	return { res.GetX(), res.GetY(), res.GetZ() };
}

XYZ_C LNLIB_XYZ_subtract(XYZ_C a, XYZ_C b)
{
	LNLib::XYZ va(a.x, a.y, a.z);
	LNLib::XYZ vb(b.x, b.y, b.z);
	LNLib::XYZ res = va - vb;
	return { res.GetX(), res.GetY(), res.GetZ() };
}

XYZ_C LNLIB_XYZ_negative(XYZ_C xyz)
{
	LNLib::XYZ va(xyz.x, xyz.y, xyz.z);
	LNLib::XYZ res = -va;
	return { res.GetX(), res.GetY(), res.GetZ() };
}

XYZ_C LNLIB_XYZ_multiply(XYZ_C xyz, double scalar)
{
	LNLib::XYZ va(xyz.x, xyz.y, xyz.z);
	LNLib::XYZ res = va * scalar;
	return { res.GetX(), res.GetY(), res.GetZ() };
}

XYZ_C LNLIB_XYZ_divide(XYZ_C xyz, double scalar)
{
	LNLib::XYZ va(xyz.x, xyz.y, xyz.z);
	LNLib::XYZ res = va / scalar;
	return { res.GetX(), res.GetY(), res.GetZ() };
}

XYZ_C LNLIB_XYZ_create_random_orthogonal(XYZ_C xyz)
{
	LNLib::XYZ vv(xyz.x, xyz.y, xyz.z);
	LNLib::XYZ res = LNLib::XYZ::CreateRandomOrthogonal(vv);
	return { res.GetX(), res.GetY(), res.GetZ() };
}

double LNLIB_XYZ_length(XYZ_C xyz)
{
	LNLib::XYZ vv(xyz.x, xyz.y, xyz.z);
	return vv.Length();
}

double LNLIB_XYZ_sqr_length(XYZ_C xyz)
{
	LNLib::XYZ vv(xyz.x, xyz.y, xyz.z);
	return vv.SqrLength();
}

double LNLIB_XYZ_angle_to(XYZ_C a, XYZ_C b, double epsilon)
{
	LNLib::XYZ va(a.x, a.y, a.z);
	LNLib::XYZ vb(b.x, b.y, b.z);
	return va.AngleTo(vb);
}

int LNLIB_XYZ_is_zero(XYZ_C xyz, double epsilon)
{
	LNLib::XYZ vv(xyz.x, xyz.y, xyz.z);
	return vv.IsZero(epsilon) ? 1 : 0;
}

int LNLIB_XYZ_is_unit(XYZ_C xyz, double epsilon)
{
	LNLib::XYZ vv(xyz.x, xyz.y, xyz.z);
	return vv.IsUnit(epsilon) ? 1 : 0;
}

int LNLIB_XYZ_is_almost_equal(XYZ_C a, XYZ_C b, double epsilon)
{
	return
		(LNLib::MathUtils::IsAlmostEqualTo(a.x, b.x, epsilon) &&
			LNLib::MathUtils::IsAlmostEqualTo(a.y, b.y, epsilon) &&
			LNLib::MathUtils::IsAlmostEqualTo(a.z, b.z, epsilon)) ? 1 : 0;
}

XYZ_C LNLIB_XYZ_normalize(XYZ_C xyz)
{
	LNLib::XYZ vv(xyz.x, xyz.y, xyz.z);
	LNLib::XYZ res = vv.Normalize();
	return { res.GetX(), res.GetY(), res.GetZ() };
}

double LNLIB_XYZ_dot_product(XYZ_C a, XYZ_C b)
{
	LNLib::XYZ va(a.x, a.y, a.z);
	LNLib::XYZ vb(b.x, b.y, b.z);
	return va.DotProduct(vb);
}

XYZ_C LNLIB_XYZ_cross_product(XYZ_C a, XYZ_C b)
{
	LNLib::XYZ va(a.x, a.y, a.z);
	LNLib::XYZ vb(b.x, b.y, b.z);
	LNLib::XYZ res = va.CrossProduct(vb);
	return { res.GetX(), res.GetY(), res.GetZ() };
}

double LNLIB_XYZ_distance(XYZ_C a, XYZ_C b)
{
	LNLib::XYZ va(a.x, a.y, a.z);
	LNLib::XYZ vb(b.x, b.y, b.z);
	return va.Distance(vb);
}

