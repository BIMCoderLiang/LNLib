/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "BezierCurve_CAPI.h"
#include "BezierCurve.h"
#include "LNObject.h"
#include <vector>

extern "C" {

	static LNLib::XYZ ToXYZ(XYZ_C c) { return LNLib::XYZ(c.x, c.y, c.z); }
	static XYZ_C FromXYZ(const LNLib::XYZ& v) { return { v.GetX(), v.GetY(), v.GetZ() }; }

	static LNLib::XYZW ToXYZW(XYZW_C c) { return LNLib::XYZW(c.wx, c.wy, c.wz, c.w); }
	static XYZW_C FromXYZW(const LNLib::XYZW& v) { return { v.GetWX(), v.GetWY(), v.GetWZ(), v.GetW() }; }

	LNLIB_EXPORT XYZ_C bezier_curve_get_point_by_bernstein(
		int degree,
		const XYZ_C* control_points,
		int cp_count,
		double paramT)
	{
		std::vector<LNLib::XYZ> cps;
		for (int i = 0; i < cp_count; ++i) cps.push_back(ToXYZ(control_points[i]));

		LNLib::LN_BezierCurve<LNLib::XYZ> curve{ degree, cps };
		LNLib::XYZ pt = LNLib::BezierCurve::GetPointOnCurveByBernstein(curve, paramT);
		return FromXYZ(pt);
	}

	LNLIB_EXPORT XYZW_C bezier_curve_get_point_by_bernstein_rational(
		int degree,
		const XYZW_C* control_points,
		int cp_count,
		double paramT)
	{
		std::vector<LNLib::XYZW> cps;
		for (int i = 0; i < cp_count; ++i) cps.push_back(ToXYZW(control_points[i]));

		LNLib::LN_BezierCurve<LNLib::XYZW> curve{ degree, cps };
		LNLib::XYZW pt = LNLib::BezierCurve::GetPointOnCurveByBernstein(curve, paramT);
		return FromXYZW(pt);
	}

	LNLIB_EXPORT XYZ_C bezier_curve_get_point_by_de_casteljau(
		int degree,
		const XYZ_C* control_points,
		int cp_count,
		double paramT)
	{
		std::vector<LNLib::XYZ> cps;
		for (int i = 0; i < cp_count; ++i) cps.push_back(ToXYZ(control_points[i]));

		LNLib::LN_BezierCurve<LNLib::XYZ> curve{ degree, cps };
		LNLib::XYZ pt = LNLib::BezierCurve::GetPointOnCurveByDeCasteljau(curve, paramT);
		return FromXYZ(pt);
	}

}