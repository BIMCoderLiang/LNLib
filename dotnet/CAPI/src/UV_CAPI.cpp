/*
 * Author:
 * 2025/11/25 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "UV_CAPI.h"
#include "UV.h"
#include "MathUtils.h"


UV_C LNLIB_UV_create(double u, double v)
{
	return { u, v };
}

double LNLIB_UV_get_u(UV_C uv)
{
	return uv.u;
}

double LNLIB_UV_get_v(UV_C uv)
{
	return uv.v;
}

UV_C LNLIB_UV_set_u(UV_C uv, double u)
{
	return { u, uv.v };
}

UV_C LNLIB_UV_set_v(UV_C uv, double v)
{
	return { uv.u, v };
}

UV_C LNLIB_UV_add(UV_C a, UV_C b)
{
	LNLib::UV va(a.u, a.v);
	LNLib::UV vb(b.u, b.v);
	LNLib::UV result = va + vb;
	return { result.GetU(), result.GetV() };
}

UV_C LNLIB_UV_subtract(UV_C a, UV_C b)
{
	LNLib::UV va(a.u, a.v);
	LNLib::UV vb(b.u, b.v);
	LNLib::UV result = va - vb;
	return { result.GetU(), result.GetV() };
}

UV_C LNLIB_UV_negative(UV_C uv)
{
	LNLib::UV v(uv.u, uv.v);
	LNLib::UV result = -v;
	return { result.GetU(), result.GetV() };
}

UV_C LNLIB_UV_normalize(UV_C uv)
{
	LNLib::UV v(uv.u, uv.v);
	LNLib::UV result = v.Normalize();
	return { result.GetU(), result.GetV() };
}

UV_C LNLIB_UV_scale(UV_C uv, double factor)
{
	LNLib::UV v(uv.u, uv.v);
	v *= factor;
	return { v.GetU(), v.GetV() };
}

UV_C LNLIB_UV_divide(UV_C uv, double divisor)
{
	if (divisor == 0.0) return uv;
	LNLib::UV v(uv.u, uv.v);
	v /= divisor;
	return { v.GetU(), v.GetV() };
}

double LNLIB_UV_length(UV_C uv)
{
	LNLib::UV v(uv.u, uv.v);
	return v.Length();
}

double LNLIB_UV_sqr_length(UV_C uv)
{
	LNLib::UV v(uv.u, uv.v);
	return v.SqrLength();
}

double LNLIB_UV_distance(UV_C a, UV_C b)
{
	LNLib::UV va(a.u, a.v);
	LNLib::UV vb(b.u, b.v);
	return va.Distance(vb);
}

int LNLIB_UV_is_zero(UV_C uv, double epsilon)
{
	LNLib::UV v(uv.u, uv.v);
	return v.IsZero(epsilon) ? 1 : 0;
}

int LNLIB_UV_is_unit(UV_C uv, double epsilon)
{
	LNLib::UV v(uv.u, uv.v);
	return v.IsUnit(epsilon) ? 1 : 0;
}

int LNLIB_UV_is_almost_equal(UV_C a, UV_C b, double epsilon)
{
	return (LNLib::MathUtils::IsAlmostEqualTo(a.u, b.u, epsilon) &&
		LNLib::MathUtils::IsAlmostEqualTo(a.v, b.v, epsilon)) ? 1 : 0;
}

double LNLIB_UV_dot_product(UV_C a, UV_C b)
{
	LNLib::UV va(a.u, a.v);
	LNLib::UV vb(b.u, b.v);
	return va.DotProduct(vb);
}

double LNLIB_UV_cross_product(UV_C a, UV_C b)
{
	LNLib::UV va(a.u, a.v);
	LNLib::UV vb(b.u, b.v);
	return va.CrossProduct(vb);
}

