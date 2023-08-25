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

std::vector<std::vector<std::vector<std::vector<LNLib::XYZ>>>> LNLib::BsplineSurface::ComputeControlPointsOfDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, int degreeU, int degreeV, UV uv, int minSpanIndexU, int maxSpanIndexU, int minSpanIndexV, int maxSpanIndexV, int derivative)
{
	std::vector<std::vector<std::vector<std::vector<LNLib::XYZ>>>> controlPointsOfDerivative(derivative + 1);
	for (int k = 0; k <= derivative; k++)
	{
		controlPointsOfDerivative[k].resize(derivative + 1);
		for (int l = 0; l <= derivative; l++)
		{
			controlPointsOfDerivative[k][l].resize(maxSpanIndexU - derivative - minSpanIndexU + 1);
			for (int i = 0; i <= maxSpanIndexU - derivative - minSpanIndexU; i++)
			{
				controlPointsOfDerivative[k][l][i].resize(maxSpanIndexV - derivative - minSpanIndexV + 1);
			}
		}
	}

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
			int dd = std::min(static_cast<int>(derivative - k), dv);
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
