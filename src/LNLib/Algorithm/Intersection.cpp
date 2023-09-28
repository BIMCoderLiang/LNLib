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
#include "LNLibExceptions.h"

using namespace LNLib;

CurveCurveIntersectionType Intersection::ComputeRays(const XYZ& point0, const XYZ& vector0, const XYZ& point1, const XYZ& vector1, double& param0, double& param1, XYZ& intersectPoint)
{
	VALIDATE_ARGUMENT(!vector0.IsZero(), "vector0", "Vector0 must not be zero vector.");
	VALIDATE_ARGUMENT(!vector1.IsZero(), "vector1", "Vector1 must not be zero vector.")

	XYZ v0 = const_cast<XYZ&>(vector0);
	XYZ v1 = const_cast<XYZ&>(vector1);

	XYZ cross = v0.CrossProduct(v1);
	XYZ diff = point1 - point0;
	XYZ coinCross = diff.CrossProduct(v1);

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

	double crossLength = cross.Length();
	double squareLength = crossLength * crossLength;

	XYZ pd1Cross = diff.CrossProduct(vector1);
	double pd1Dot = pd1Cross.DotProduct(cross);
	param0 = pd1Dot / squareLength;

	XYZ pd2Cross = diff.CrossProduct(vector0);
	double pd2Dot = pd2Cross.DotProduct(cross);
	param1 = pd2Dot / squareLength;

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
