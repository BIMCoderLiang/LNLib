/**
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "BsplineCurve.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "ValidationUtils.h"
#include <math.h>
#include <algorithm>


int LNLib::BsplineCurve::GetContinuity(unsigned int degree, const std::vector<double>& knotVector, double knot)
{
	int multi = Polynomials::GetKnotMultiplicity(knot, knotVector);
	return degree - multi;
}

void LNLib::BsplineCurve::ComputeControlPointsOfDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, unsigned int derivative, unsigned int derMin, unsigned int derMax, std::vector<std::vector<XYZ>>& controlPointsOfDerivative)
{
	controlPointsOfDerivative.resize(derivative + 1);
	for (int i = 0; i <= static_cast<int>(derivative); i++)
	{
		controlPointsOfDerivative[i].resize(derMax - derivative - derMin + 1);
	}

	int range = derMax - derMin;
	for (int i = 0; i <= range; i++)
	{
		controlPointsOfDerivative[0][i] = controlPoints[derMin + i];
	}
	for (int k = 1; k <= static_cast<int>(derivative); k++)
	{
		int temp = degree - k + 1;
		for (int i = 0; i <= range - k; i++)
		{
			controlPointsOfDerivative[k][i] = temp * (controlPointsOfDerivative[k - 1][i + 1] - controlPointsOfDerivative[k - 1][i]) / (knotVector[derMin + i + degree + 1] - knotVector[derMin + i + k]);
		}
	}
}
