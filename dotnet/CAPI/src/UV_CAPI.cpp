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

extern "C" {

	LNLIB_EXPORT UV_C uv_create(double u, double v)
	{
		return { u, v };
	}

	LNLIB_EXPORT double uv_get_u(UV_C uv)
	{
		return uv.u;
	}

	LNLIB_EXPORT double uv_get_v(UV_C uv)
	{
		return uv.v;
	}

	LNLIB_EXPORT UV_C uv_add(UV_C a, UV_C b)
	{
		LNLib::UV va(a.u, a.v);
		LNLib::UV vb(b.u, b.v);
		LNLib::UV result = va + vb;
		return { result.GetU(), result.GetV() };
	}

	LNLIB_EXPORT UV_C uv_subtract(UV_C a, UV_C b)
	{
		LNLib::UV va(a.u, a.v);
		LNLib::UV vb(b.u, b.v);
		LNLib::UV result = va - vb;
		return { result.GetU(), result.GetV() };
	}

	LNLIB_EXPORT UV_C uv_negative(UV_C uv)
	{
		LNLib::UV v(uv.u, uv.v);
		LNLib::UV result = -v;
		return { result.GetU(), result.GetV() };
	}

	LNLIB_EXPORT UV_C uv_normalize(UV_C uv)
	{
		LNLib::UV v(uv.u, uv.v);
		LNLib::UV result = v.Normalize();
		return { result.GetU(), result.GetV() };
	}

	LNLIB_EXPORT UV_C uv_scale(UV_C uv, double factor)
	{
		LNLib::UV v(uv.u, uv.v);
		v *= factor;
		return { v.GetU(), v.GetV() };
	}

	LNLIB_EXPORT UV_C uv_divide(UV_C uv, double divisor)
	{
		if (divisor == 0.0) return uv;
		LNLib::UV v(uv.u, uv.v);
		v /= divisor;
		return { v.GetU(), v.GetV() };
	}

	LNLIB_EXPORT double uv_length(UV_C uv)
	{
		LNLib::UV v(uv.u, uv.v);
		return v.Length();
	}

	LNLIB_EXPORT double uv_sqr_length(UV_C uv)
	{
		LNLib::UV v(uv.u, uv.v);
		return v.SqrLength();
	}

	LNLIB_EXPORT double uv_distance(UV_C a, UV_C b)
	{
		LNLib::UV va(a.u, a.v);
		LNLib::UV vb(b.u, b.v);
		return va.Distance(vb);
	}

	LNLIB_EXPORT int uv_is_zero(UV_C uv, double epsilon)
	{
		LNLib::UV v(uv.u, uv.v);
		return v.IsZero(epsilon) ? 1 : 0;
	}

	LNLIB_EXPORT int uv_is_unit(UV_C uv, double epsilon)
	{
		LNLib::UV v(uv.u, uv.v);
		return v.IsUnit(epsilon) ? 1 : 0;
	}

	LNLIB_EXPORT int uv_is_almost_equal(UV_C a, UV_C b, double epsilon)
	{
		return (LNLib::MathUtils::IsAlmostEqualTo(a.u, b.u, epsilon) &&
			LNLib::MathUtils::IsAlmostEqualTo(a.v, b.v, epsilon)) ? 1 : 0;
	}

	LNLIB_EXPORT double uv_dot(UV_C a, UV_C b)
	{
		LNLib::UV va(a.u, a.v);
		LNLib::UV vb(b.u, b.v);
		return va.DotProduct(vb);
	}

	LNLIB_EXPORT double uv_cross(UV_C a, UV_C b)
	{
		LNLib::UV va(a.u, a.v);
		LNLib::UV vb(b.u, b.v);
		return va.CrossProduct(vb);
	}

}
