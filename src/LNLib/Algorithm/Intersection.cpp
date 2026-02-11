/*
 * Author:
 * 2023/06/23 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Intersection.h"
#include "XYZ.h"
#include "MathUtils.h"
#include "LNLibExceptions.h"

using namespace LNLib;

CurveCurveIntersectionType Intersection::ComputeRays(
    const XYZ& point0, const XYZ& vector0,
    const XYZ& point1, const XYZ& vector1,
    double& param0, double& param1, XYZ& intersectPoint)
{
    if (vector0.IsAlmostEqualTo(vector1))
    {
        if (point0.IsAlmostEqualTo(point1))
        {
            intersectPoint = point0;
            param0 = param1 = 0;
            return CurveCurveIntersectionType::Intersecting;
        }
    }

    VALIDATE_ARGUMENT(!vector0.IsZero(), "vector0", "Vector0 must not be zero.");
    VALIDATE_ARGUMENT(!vector1.IsZero(), "vector1", "Vector1 must not be zero.");

    XYZ diff = point0 - point1;

    double a = vector0.DotProduct(vector0);
    double b = vector0.DotProduct(vector1);
    double c = vector1.DotProduct(vector1);
    double d = vector0.DotProduct(diff);
    double e = vector1.DotProduct(diff);

    double denom = a * c - b * b;

    if (MathUtils::IsAlmostEqualTo(std::abs(denom), 0.0)) {
        XYZ w = point1 - point0;
        if (MathUtils::IsAlmostEqualTo(w.CrossProduct(vector0).Length(), 0.0)) {
            param0 = param1 = 0.0;
            intersectPoint = point0;
            return CurveCurveIntersectionType::Coincident;
        }
        else {
            return CurveCurveIntersectionType::Parallel;
        }
    }

    param0 = (b * e - c * d) / denom;
    param1 = (a * e - b * d) / denom;

    XYZ p0 = point0 + vector0 * param0;
    XYZ p1 = point1 + vector1 * param1;

    if (MathUtils::IsAlmostEqualTo(p0.Distance(p1),0.0)) {
        intersectPoint = (p0 + p1) * 0.5;
        return CurveCurveIntersectionType::Intersecting;
    }
    else {
        intersectPoint = p0; 
        return CurveCurveIntersectionType::Skew;
    }
}

LinePlaneIntersectionType LNLib::Intersection::ComputeLineAndPlane(
	const XYZ& normal,
	const XYZ& pointOnPlane,
	const XYZ& pointOnLine,
	const XYZ& lineDirection,
	XYZ& intersectPoint)
{
	double denom = lineDirection.DotProduct(normal);

	XYZ P0Q = pointOnLine - pointOnPlane;
	double numer = -P0Q.DotProduct(normal);

	if (MathUtils::IsAlmostEqualTo(std::abs(denom),0.0)) {
		if (MathUtils::IsAlmostEqualTo(std::abs(numer), 0.0)) {
			intersectPoint = pointOnLine;
			return LinePlaneIntersectionType::On;
		}
		else {
			return LinePlaneIntersectionType::Parallel;
		}
	}

	double t = numer / denom;
	intersectPoint = pointOnLine + t * lineDirection;
	return LinePlaneIntersectionType::Intersecting;
}
