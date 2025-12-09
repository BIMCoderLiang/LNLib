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

XYZ_C LNLIB_BEZIERCUR_get_point_on_curve_by_bernstein(
	int degree,
	const XYZ_C* control_points,
	int control_points_count,
	double param_t)
{
	std::vector<LNLib::XYZ> cps;
	for (int i = 0; i < control_points_count; ++i) cps.push_back(ToXYZ(control_points[i]));

	LNLib::LN_BezierCurve<LNLib::XYZ> curve{ degree, cps };
	LNLib::XYZ pt = LNLib::BezierCurve::GetPointOnCurveByBernstein(curve, param_t);
	return FromXYZ(pt);
}

XYZW_C LNLIB_BEZIERCUR_get_rational_point_on_curve_by_bernstein(
	int degree,
	const XYZW_C* control_points,
	int control_points_count,
	double param_t)
{
	std::vector<LNLib::XYZW> cps;
	for (int i = 0; i < control_points_count; ++i) cps.push_back(ToXYZW(control_points[i]));

	LNLib::LN_BezierCurve<LNLib::XYZW> curve{ degree, cps };
	LNLib::XYZW pt = LNLib::BezierCurve::GetPointOnCurveByBernstein(curve, param_t);
	return FromXYZW(pt);
}

XYZ_C LNLIB_BEZIERCUR_get_point_on_curve_by_deCasteljau(
	int degree,
	const XYZ_C* control_points,
	int control_points_count,
	double param_t)
{
	std::vector<LNLib::XYZ> cps;
	for (int i = 0; i < control_points_count; ++i) cps.push_back(ToXYZ(control_points[i]));

	LNLib::LN_BezierCurve<LNLib::XYZ> curve{ degree, cps };
	LNLib::XYZ pt = LNLib::BezierCurve::GetPointOnCurveByDeCasteljau(curve, param_t);
	return FromXYZ(pt);
}

XYZW_C LNLIB_BEZIERCUR_get_rational_point_on_curve_by_deCasteljau(
	int degree,
	const XYZW_C* control_points,
	int control_points_count,
	double param_t)
{
	std::vector<LNLib::XYZW> cps;
	for (int i = 0; i < control_points_count; ++i) cps.push_back(ToXYZW(control_points[i]));

	LNLib::LN_BezierCurve<LNLib::XYZW> curve{ degree, cps };
	LNLib::XYZW pt = LNLib::BezierCurve::GetPointOnCurveByDeCasteljau(curve, param_t);
	return FromXYZW(pt);
}