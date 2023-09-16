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

double LNLib::ValidationUtils::ComputeMaxErrorOfBezierReduction(int degree, const std::vector<XYZW>& currentControlPoints, std::vector<XYZW>& reductedControlPoints)
{
	int r = floor((degree - 1)/2);
	
	reductedControlPoints[0] = currentControlPoints[0];
	if (degree % 2 == 0)
	{
		for (int i = 1; i <= r; i++)
		{
			double alpha = (double)i / (double)degree;
			reductedControlPoints[i] = (currentControlPoints[i] - alpha * reductedControlPoints[i - 1]) / (1 - alpha);
		}
		for (int i = degree - 2; i <= r + 1; i++)
		{
			double alpha = (double)(i + 1) / (double)degree;
			reductedControlPoints[i] = (currentControlPoints[i + 1] - (1 - alpha) * reductedControlPoints[i + 1]) / alpha;
		}

		reductedControlPoints[degree - 1] = currentControlPoints[degree];
		return currentControlPoints[r + 1].Distance(0.5 * (reductedControlPoints[r] + reductedControlPoints[r + 1]));
	}
	else
	{
		for (int i = 1; i <= r - 1; i++)
		{
			double alpha = (double)i / (double)degree;
			reductedControlPoints[i] = (currentControlPoints[i] - alpha * reductedControlPoints[i - 1]) / (1 - alpha);
		}
		for (int i = degree - 2; i <= r; i++)
		{
			double alpha = (double)(i + 1) / (double)degree;
			reductedControlPoints[i] = (currentControlPoints[i + 1] - (1 - alpha) * reductedControlPoints[i + 1]) / alpha;
		}

		double alpha = (double)r / (double)degree;
		XYZW PLr = (currentControlPoints[r] - alpha * reductedControlPoints[r - 1]) / (1 - alpha);
		alpha = (double)(r + 1) / (double)degree;
		XYZW PRr = (currentControlPoints[r + 1] - (1 - alpha) * reductedControlPoints[r + 1]) / alpha;
		reductedControlPoints[r] = 0.5 * (PLr + PRr);
		reductedControlPoints[degree - 1] = currentControlPoints[degree];
		return std::abs((PLr.ToXYZ(true) + PRr.ToXYZ(true)).Length());
	}
}



