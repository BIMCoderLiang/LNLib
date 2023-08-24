/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "BsplineCurve.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "ValidationUtils.h"
#include <math.h>
#include <algorithm>


int LNLib::BsplineCurve::GetContinuity(int degree, const std::vector<double>& knotVector, double knot)
{
	int multi = Polynomials::GetKnotMultiplicity(knotVector, knot);
	return degree - multi;
}

std::vector<std::vector<LNLib::XYZ>> LNLib::BsplineCurve::ComputeControlPointsOfDerivatives(int degree, int derivative, int minSpanIndex, int maxSpanIndex, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	int range = maxSpanIndex - minSpanIndex;
	std::vector<std::vector<LNLib::XYZ>> PK(derivative + 1, std::vector<XYZ>(range + 1));

	for (int i = 0; i <= range; i++)
	{
		PK[0][i] = controlPoints[minSpanIndex + i];
	}
	for (int k = 1; k <= derivative; k++)
	{
		int temp = degree - k + 1;
		for (int i = 0; i <= range - k; i++)
		{
			PK[k][i] = temp * (PK[k - 1][i + 1] - PK[k - 1][i]) / (knotVector[minSpanIndex + i + degree + 1] - knotVector[minSpanIndex + i + k]);
		}
	}
	return PK;
}
