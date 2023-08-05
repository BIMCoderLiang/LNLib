/*
 * Author:
 * 2023/06/23 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "Intersection.h"
#include "XYZ.h"
#include "MathUtils.h"

using namespace LNLib;

CurveCurveIntersectionType Intersection::ComputeRays(const XYZ& point0, const XYZ& vector0, const XYZ& point1, const XYZ& vector1, double& param0, double& param1, XYZ& intersectPoint)
{
	XYZ v0 = const_cast<XYZ&>(vector0).Normalize();
	XYZ v1 = const_cast<XYZ&>(vector1).Normalize();

	XYZ cross = v0.CrossProduct(v1);

	XYZ p0Temp = point0;
	XYZ p1Temp = point1;

	XYZ pDiff = p1Temp - p0Temp;
	XYZ pDiffNormal = (p1Temp - p0Temp).Normalize();
	XYZ coinCross = pDiffNormal.CrossProduct(v1);

	if (cross.IsZero())
	{
		if (coinCross.IsZero())
		{
			return CurveCurveIntersectionType::Coincident;
		}
		else
		{
			return CurveCurveIntersectionType::Parallel;
		}
	}

	double dMagn = cross.Length();
	double dMagnSquare = dMagn * dMagn;

	XYZ pd0Cross = pDiff.Normalize().CrossProduct(vector1);
	double d0 = pd0Cross.DotProduct(cross);
	param0 = d0 / dMagnSquare;

	XYZ pd1Cross = pDiff.Normalize().CrossProduct(vector0);
	double d1 = pd1Cross.DotProduct(cross);
	param1 = d1 / dMagnSquare;

	XYZ rayP0 = point0 + vector0 * param0;
	XYZ rayP1 = point1 + vector1 * param1;

	if (rayP0.IsAlmostEqualTo(rayP1))
	{
		intersectPoint = rayP0;
		return CurveCurveIntersectionType::Intersecting;
	}
	return CurveCurveIntersectionType::Skew;
}

LinePlaneIntersectionType LNLib::Intersection::ComputeLineAndPlane(const XYZ& normal, const XYZ& pointOnPlane, const XYZ& pointOnLine, const XYZ& lineDirection, XYZ& intersectPoint)
{
	double d = (pointOnLine - pointOnPlane).DotProduct(normal);
	double angle = normal.AngleTo(lineDirection);
	if (MathUtils::IsAlmostEqualTo(d, 0.0))
	{
		return LinePlaneIntersectionType::On;
	}
	if (MathUtils::IsAlmostEqualTo(angle, Constants::Pi / 2))
	{
		return LinePlaneIntersectionType::Parallel;
	}
	d = d / lineDirection.DotProduct(normal);
	intersectPoint = d * const_cast<XYZ&>(normal).Normalize() + pointOnLine;
	return LinePlaneIntersectionType::Intersecting;
}
