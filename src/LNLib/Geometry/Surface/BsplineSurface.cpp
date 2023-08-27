/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */


#include "BsplineSurface.h"
#include "BsplineCurve.h"
#include "Polynomials.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "ValidationUtils.h"
#include <algorithm>

std::vector<std::vector<std::vector<std::vector<LNLib::XYZ>>>> LNLib::BsplineSurface::ComputeControlPointsOfDerivatives(int degreeU, int degreeV, int derivative, int minSpanIndexU, int maxSpanIndexU, int minSpanIndexV, int maxSpanIndexV, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, UV uv, const std::vector<std::vector<XYZ>>& controlPoints)
{
	VALIDATE_ARGUMENT(degreeU > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(degreeV > 0, "degreeU", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(derivative >= 0, "derivative", "derivative must greater than or equals zero.");
	VALIDATE_ARGUMENT_RANGE(minSpanIndexU, 0, maxSpanIndexU);
	VALIDATE_ARGUMENT_RANGE(minSpanIndexV, 0, maxSpanIndexV);
	VALIDATE_ARGUMENT(knotVectorU.size() > 0, "knotVectorU", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorU), "knotVectorU", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetU(), knotVectorU[0], knotVectorU[knotVectorU.size() - 1]);
	VALIDATE_ARGUMENT(knotVectorV.size() > 0, "knotVectorV", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVectorV), "knotVectorV", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(uv.GetV(), knotVectorV[0], knotVectorV[knotVectorV.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeU, knotVectorU.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidBspline(degreeV, knotVectorV.size(), controlPoints[0].size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<std::vector<std::vector<LNLib::XYZ>>>> controlPointsOfDerivative(derivative + 1, 
		std::vector<std::vector<std::vector<LNLib::XYZ>>>(maxSpanIndexU - minSpanIndexU + 1, 
			std::vector<std::vector<LNLib::XYZ>>(maxSpanIndexV - minSpanIndexV + 1)));

	int du = std::min(derivative, degreeU);
	int dv = std::min(derivative, degreeV);
	int rangeU = maxSpanIndexU - minSpanIndexU;
	int rangeV = maxSpanIndexV - minSpanIndexV;

	for (int j = minSpanIndexV; j <= maxSpanIndexV; j++)
	{
		std::vector<XYZ> points;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			points.emplace_back(controlPoints[i][j]);
		}

		std::vector<std::vector<XYZ>> temp = BsplineCurve::ComputeControlPointsOfDerivatives(degreeU, du, minSpanIndexU, maxSpanIndexU, knotVectorU, points);
		for (int k = 0; k <= du; k++)
		{
			for (int i = 0; i <= rangeU - k; i++)
			{
				controlPointsOfDerivative[k][0][i][j - minSpanIndexV] = temp[k][i];
			}
		}
	}
	for (int k = 0; k < du; k++)
	{
		for (int i = 0; i <= rangeU - k; i++)
		{
			std::vector<XYZ> points = controlPointsOfDerivative[k][0][i];
			int dd = std::min(derivative - k, dv);
			std::vector<std::vector<XYZ>> temp = BsplineCurve::ComputeControlPointsOfDerivatives(degreeV, dd, 0, rangeV, knotVectorV, points);
			for (int l = 1; l <= dd; l++)
			{
				for (int j = 0; j < rangeV - l; j++)
				{
					controlPointsOfDerivative[k][l][i][j] = temp[l][j];
				}
			}
		}
	}
	return controlPointsOfDerivative;
}
