/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "ValidationUtils.h"
#include "XYZW.h"
#include "Constants.h"
#include "LNLibExceptions.h"
#include <algorithm>


namespace LNLib
{
	static double GetCoefficient(int index, int degree)
	{
		return index / degree;
	}
}

bool LNLib::ValidationUtils::IsValidBezier(int degree, int controlPointsCount)
{
	return controlPointsCount == degree + 1;
}

bool LNLib::ValidationUtils::IsValidKnotVector(const std::vector<double>& knotVector)
{
	return std::is_sorted(knotVector.begin(), knotVector.end());
}

bool LNLib::ValidationUtils::IsValidBspline(int degree, int knotVectorCount, int controlPointsCount)
{
	return (knotVectorCount - 1) == (controlPointsCount - 1) + degree + 1;
}

bool LNLib::ValidationUtils::IsValidNurbs(int degree, int knotVectorCount, int weightedControlPointsCount)
{
	return (knotVectorCount - 1) == (weightedControlPointsCount - 1) + degree + 1;
}

bool LNLib::ValidationUtils::IsValidDegreeReduction(int degree)
{
	return degree > 1;
}

double LNLib::ValidationUtils::ComputeCurveModifyTolerance(const std::vector<XYZW>& controlPoints)
{
	double minWeight = 1.0;
	double maxDistance = 0.0;

	int size = static_cast<int>(controlPoints.size());
	for (int i = 0; i < size; i++)
	{
		XYZW temp = controlPoints[i];
		minWeight = std::min(minWeight, temp.GetW());
		maxDistance = std::max(maxDistance, temp.ToXYZ(true).Length());
	}

	return Constants::DistanceEpsilon * minWeight / (1 + std::abs(maxDistance));
}

double LNLib::ValidationUtils::ComputeMaxErrorOfBezierReduction(int degree, const std::vector<XYZW>& currentControlPoints, const std::vector<XYZW>& reductedControlPoints)
{
	int r = (degree - 1)/2;

	if (degree % 2 == 0)
	{
		XYZW cwp = currentControlPoints[r];
		XYZW rwp = reductedControlPoints[r-1];
		XYZ Pr = (cwp.ToXYZ(true) - GetCoefficient(r, degree) * rwp.ToXYZ(true)) / (1 - GetCoefficient(r, degree));

		XYZW cwp1 = currentControlPoints[r + 2];
		XYZW rwp1 = reductedControlPoints[r + 2];
		XYZ Pr1 = (cwp1.ToXYZ(true) - (1 - GetCoefficient(r + 2, degree)) * rwp1.ToXYZ(true)) / GetCoefficient(r + 2,degree);

		XYZW cp1 = currentControlPoints[r + 1];
		XYZ pr = 0.5 * (Pr + Pr1);

		return std::abs((cp1.ToXYZ(true)-pr).Length());
	}
	else
	{
		XYZW cwp = currentControlPoints[r];
		XYZW rwp = reductedControlPoints[r - 1];
		XYZ PLr = (cwp.ToXYZ(true) - GetCoefficient(r, degree) * rwp.ToXYZ(true)) / (1 - GetCoefficient(r, degree));

		XYZW cwp1 = currentControlPoints[r + 1];
		XYZW rwp1 = reductedControlPoints[r + 1];
		XYZ PRr = (cwp1.ToXYZ(true) - (1 - GetCoefficient(r + 1, degree)) * rwp1.ToXYZ(true)) / GetCoefficient(r + 1, degree);

		return std::abs((PLr+PRr).Length());
	}
}



