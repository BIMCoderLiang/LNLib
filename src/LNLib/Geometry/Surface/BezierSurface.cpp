/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "BezierSurface.h"
#include "BezierCurve.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "LNLibExceptions.h"

using namespace LNLib;

XYZ BezierSurface::GetPointOnSurfaceByDeCasteljau(int degreeU, int degreeV, const std::vector<std::vector<XYZ>>& controlPoints, UV uv)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(degreeU + 1 == static_cast<int>(controlPoints.size()), "controlPoints", "ControlPoints row size equals degreeU plus one.");
	VALIDATE_ARGUMENT(degreeV + 1 == static_cast<int>(controlPoints[0].size()), "controlPoints", "ControlPoints column size equals degreeV plus one.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), 0.0, 1.0);
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), 0.0, 1.0);

	XYZ point(0,0,0);
	std::vector<XYZ> temp(degreeU + 1);
	for (int i = 0; i <= degreeU; i++)
	{
		temp[i] = BezierCurve::GetPointOnCurveByDeCasteljau(degreeV, controlPoints[i], uv.GetV());
	}
	point = BezierCurve::GetPointOnCurveByDeCasteljau(degreeU, temp, uv.GetU());
	return point;
}

XYZW LNLib::BezierSurface::GetPointOnRationalSurfaceByDeCasteljau(int degreeU, int degreeV, const std::vector<std::vector<XYZW>>& controlPoints, UV uv)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeV", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(degreeU + 1 == static_cast<int>(controlPoints.size()), "controlPoints", "ControlPoints row size equals degreeU plus one.");
	VALIDATE_ARGUMENT(degreeV + 1 == static_cast<int>(controlPoints[0].size()), "controlPoints", "ControlPoints column size equals degreeV plus one.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), 0.0, 1.0);
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), 0.0, 1.0);

	XYZW point(0, 0, 0, 0);
	std::vector<XYZW> temp(degreeU + 1);
	for (int i = 0; i <= degreeU; i++)
	{
		temp[i] = BezierCurve::GetPointOnRationalCurveByDeCasteljau(degreeV, controlPoints[i], uv.GetV());
	}
	point = BezierCurve::GetPointOnRationalCurveByDeCasteljau(degreeU, temp, uv.GetU());
	return point;
}
