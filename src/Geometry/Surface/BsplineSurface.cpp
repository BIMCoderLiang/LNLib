
#include "BsplineSurface.h"
#include "BsplineCurve.h"
#include "Polynomials.h"
#include "UV.h"
#include "XYZ.h"
#include "XYZW.h"
#include "ValidationUtils.h"
#include <algorithm>

void LNLib::BsplineSurface::ComputeControlPointsOfDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, UV uv, unsigned int derMinU, unsigned int derMaxU, unsigned int derMinV, unsigned int derMaxV, unsigned int derivative,  std::vector<std::vector<std::vector<std::vector<XYZ>>>>& controlPointsOfDerivative)
{
	controlPointsOfDerivative.resize(derivative + 1);
	for (int k = 0; k <= static_cast<int>(derivative); k++)
	{
		controlPointsOfDerivative[k].resize(derivative + 1);
		for (int l = 0; l <= static_cast<int>(derivative); l++)
		{
			controlPointsOfDerivative[k][l].resize(derMaxU - derivative - derMinU + 1);
			for (int i = 0; i <= static_cast<int>(derMaxU - derivative - derMinU); i++)
			{
				controlPointsOfDerivative[k][l][i].resize(derMaxV - derivative - derMinV + 1);
			}
		}
	}

	int du = std::min(derivative, degreeU);
	int dv = std::min(derivative, degreeV);
	int rangeU = derMaxU - derMinU;
	int rangeV = derMaxV - derMinV;

	for (int j = derMinV; j <= static_cast<int>(derMaxV); j++)
	{
		std::vector<XYZ> points;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			points.emplace_back(controlPoints[i][j]);
		}

		std::vector<std::vector<XYZ>> temp;
		BsplineCurve::ComputeControlPointsOfDerivatives(degreeU, knotVectorU, points, du, derMinU, derMaxU, temp);
		for (int k = 0; k <= du; k++)
		{
			for (int i = 0; i <= rangeU - k; i++)
			{
				controlPointsOfDerivative[k][0][i][j - derMinV] = temp[k][i];
			}
		}
	}
	for (int k = 0; k < du; k++)
	{
		for (int i = 0; i <= rangeU - k; i++)
		{
			std::vector<XYZ> points = controlPointsOfDerivative[k][0][i];
			int dd = std::min(static_cast<int>(derivative - k), dv);
			std::vector<std::vector<XYZ>> temp;
			BsplineCurve::ComputeControlPointsOfDerivatives(degreeV, knotVectorV, points, dd,0,rangeV,temp);
			for (int l = 1; l <= dd; l++)
			{
				for (int j = 0; j < rangeV - l; j++)
				{
					controlPointsOfDerivative[k][l][i][j] = temp[l][j];
				}
			}
		}
	}
}
