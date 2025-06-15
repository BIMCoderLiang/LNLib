/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
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

	int size = controlPoints.size();
	for (int i = 0; i < size; i++)
	{
		XYZW temp = controlPoints[i];
		minWeight = std::min(minWeight, temp.GetW());
		maxDistance = std::max(maxDistance, temp.ToXYZ(true).Length());
	}

	return Constants::DistanceEpsilon * minWeight / (1 + std::abs(maxDistance));
}



