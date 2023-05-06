
#include "BsplineSurface.h"
#include "BsplineCurve.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "ValidationUtils.h"
#include <algorithm>

void LNLib::BsplineSurface::GetPointOnSurface(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, XYZ& point)
{
	unsigned int n = static_cast<int>(controlPoints.size() - 1);
	unsigned int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, paramU, knotVectorU);
	std::vector<double> basisFunctionsU;
	basisFunctionsU.resize(degreeU + 1);
	Polynomials::BasisFunctions(uSpanIndex, degreeU, paramU, knotVectorU, basisFunctionsU);

	unsigned int m = static_cast<int>(controlPoints[0].size() - 1);
	unsigned int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, paramV, knotVectorV);
	std::vector<double> basisFunctionsV;
	basisFunctionsV.resize(degreeV + 1);
	Polynomials::BasisFunctions(vSpanIndex, degreeV, paramV, knotVectorV, basisFunctionsV);

	unsigned int uind = uSpanIndex - degreeU;
	for (unsigned int l = 0; l <= degreeV; l++)
	{
		XYZ temp = XYZ(0, 0, 0);
		int vind = vSpanIndex - degreeV + 1;
		for (unsigned int k = 0; k <= degreeU; k++)
		{
			temp += basisFunctionsU[k] * controlPoints[uind + k][vind];
		}
		point += basisFunctionsV[l] * temp;
	}
}

void LNLib::BsplineSurface::ComputeDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives)
{
	unsigned int du = std::min(derivative, degreeU);
	for (unsigned int k = degreeU + 1; k <= derivative; k++)
	{
		for (unsigned int l = 0; l <= derivative - k; l++)
		{
			derivatives[k][l] = XYZ();
		}
	}

	unsigned int dv = std::min(derivative, degreeV);
	for (unsigned int l = degreeV + 1; l <= derivative; l++)
	{
		for (unsigned int k = 0; k <= derivative - l; k++)
		{
			derivatives[k][l] = XYZ();
		}
	}

	unsigned int n = static_cast<int>(controlPoints.size() - 1);
	unsigned int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, paramU, knotVectorU);
	std::vector<std::vector<double>> derivativeBasisFunctionsU;
	derivativeBasisFunctionsU.resize(derivative + 1);
	for (unsigned i = 0; i <= derivative; i++)
	{
		derivativeBasisFunctionsU[i].resize(degreeU + 1);
	}
	Polynomials::BasisFunctionsDerivatives(uSpanIndex, degreeU, paramU, du, knotVectorU, derivativeBasisFunctionsU);

	unsigned int m = static_cast<int>(controlPoints[0].size() - 1);
	unsigned int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, paramV, knotVectorV);
	std::vector<std::vector<double>> derivativeBasisFunctionsV;
	derivativeBasisFunctionsV.resize(derivative + 1);
	for (unsigned i = 0; i <= derivative; i++)
	{
		derivativeBasisFunctionsV[i].resize(degreeV + 1);
	}
	Polynomials::BasisFunctionsDerivatives(vSpanIndex, degreeV, paramV, dv, knotVectorV, derivativeBasisFunctionsV);

	for (unsigned k = 0; k <= du; k++)
	{
		for (unsigned s = 0; s <= degreeV; s++)
		{
			XYZ temp = XYZ();
			for (unsigned r = 0; r <= degreeU; r++)
			{
				temp += derivativeBasisFunctionsU[k][r]* controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
			}
			unsigned int dd = std::min(derivative - k, dv);
			for (unsigned l = 0; l <= dd; l++)
			{
				derivatives[k][l] = XYZ();
				for (unsigned int s = 0; s <= degreeV; s++)
				{
					derivatives[k][l] += derivativeBasisFunctionsV[l][s] * temp;
				}
			}
		}
	}
}

void LNLib::BsplineSurface::ComputeDerivatives(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, std::vector<std::vector<XYZW>>& derivatives)
{
	unsigned int du = std::min(derivative, degreeU);
	for (unsigned int k = degreeU + 1; k <= derivative; k++)
	{
		for (unsigned int l = 0; l <= derivative - k; l++)
		{
			derivatives[k][l] = XYZW();
		}
	}

	unsigned int dv = std::min(derivative, degreeV);
	for (unsigned int l = degreeV + 1; l <= derivative; l++)
	{
		for (unsigned int k = 0; k <= derivative - l; k++)
		{
			derivatives[k][l] = XYZW();
		}
	}

	unsigned int n = static_cast<int>(controlPoints.size() - 1);
	unsigned int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, paramU, knotVectorU);
	std::vector<std::vector<double>> derivativeBasisFunctionsU;
	derivativeBasisFunctionsU.resize(derivative + 1);
	for (unsigned i = 0; i <= derivative; i++)
	{
		derivativeBasisFunctionsU[i].resize(degreeU + 1);
	}
	Polynomials::BasisFunctionsDerivatives(uSpanIndex, degreeU, paramU, du, knotVectorU, derivativeBasisFunctionsU);

	unsigned int m = static_cast<int>(controlPoints[0].size() - 1);
	unsigned int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, paramV, knotVectorV);
	std::vector<std::vector<double>> derivativeBasisFunctionsV;
	derivativeBasisFunctionsV.resize(derivative + 1);
	for (unsigned i = 0; i <= derivative; i++)
	{
		derivativeBasisFunctionsV[i].resize(degreeV + 1);
	}
	Polynomials::BasisFunctionsDerivatives(vSpanIndex, degreeV, paramV, dv, knotVectorV, derivativeBasisFunctionsV);

	for (unsigned k = 0; k <= du; k++)
	{
		for (unsigned s = 0; s <= degreeV; s++)
		{
			XYZW temp = XYZW();
			for (unsigned r = 0; r <= degreeU; r++)
			{
				temp += derivativeBasisFunctionsU[k][r] * controlPoints[uSpanIndex - degreeU + r][vSpanIndex - degreeV + s];
			}
			unsigned int dd = std::min(derivative - k, dv);
			for (unsigned l = 0; l <= dd; l++)
			{
				derivatives[k][l] = XYZW();
				for (unsigned int s = 0; s <= degreeV; s++)
				{
					derivatives[k][l] += derivativeBasisFunctionsV[l][s] * temp;
				}
			}
		}
	}
}

void LNLib::BsplineSurface::ComputeControlPointsOfDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, unsigned int minU, unsigned int maxU, unsigned int minV, unsigned int maxV, std::vector<std::vector<std::vector<std::vector<XYZ>>>>& controlPointsOfDerivative)
{
	int du = std::min(derivative, degreeU);
	int dv = std::min(derivative, degreeV);
	int rangeU = maxU - minU;
	int rangeV = maxV - minV;

	for (unsigned int j = minV; j <= maxV; j++)
	{
		std::vector<XYZ> points;
		for (int i = 0; i < controlPoints.size(); i++)
		{
			points.emplace_back(controlPoints[i][j]);
		}

		std::vector<std::vector<XYZ>> temp;
		BsplineCurve::ComputeControlPointsOfDerivatives(degreeU, knotVectorU, points, du, minU, maxU, temp);
		for (int k = 0; k <= du; k++)
		{
			for (int i = 0; i <= rangeU - k; i++)
			{
				controlPointsOfDerivative[k][0][i][j - minV] = temp[k][i];
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

void LNLib::BsplineSurface::ComputeDerivativesByAllBasisFunctions(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives)
{
	int du = std::min(derivative, degreeU);
	for (unsigned int k = degreeU + 1; k <= derivative; k++)
	{
		for (unsigned int l = 0; l <= derivative - k; l++)
		{
			derivatives[k][l] = XYZ();
		}
	}
	int dv = std::min(derivative, degreeV);
	for (unsigned int l = degreeV + 1; l <= derivative; l++)
	{
		for (unsigned int k = 0; k <= derivative - l; k++)
		{
			derivatives[k][l] = XYZ();
		}
	}

	unsigned int n = static_cast<int>(controlPoints.size() - 1);
	unsigned int uSpanIndex = Polynomials::GetKnotSpanIndex(n, degreeU, paramU, knotVectorU);
	std::vector<std::vector<double>> allBasisFunctionsU;
	for (unsigned int i = 0; i <= degreeU; i++)
	{
		std::vector<double> basisFunctions;
		Polynomials::BasisFunctions(uSpanIndex, i, paramU, knotVectorU, basisFunctions);
		allBasisFunctionsU[i] = basisFunctions;
	}

	unsigned int m = static_cast<int>(controlPoints[0].size() - 1);
	unsigned int vSpanIndex = Polynomials::GetKnotSpanIndex(m, degreeV, paramV, knotVectorV);
	std::vector<std::vector<double>> allBasisFunctionsV;
	for (unsigned int i = 0; i <= degreeV; i++)
	{
		std::vector<double> basisFunctions;
		Polynomials::BasisFunctions(vSpanIndex, i, paramV, knotVectorV, basisFunctions);
		allBasisFunctionsV[i] = basisFunctions;
	}

	std::vector<std::vector<std::vector<std::vector<XYZ>>>> controlPointsOfDerivative;
	ComputeControlPointsOfDerivatives(controlPoints, knotVectorU, knotVectorV, degreeU, degreeV, paramU, paramV, derivative, uSpanIndex - degreeU, uSpanIndex, vSpanIndex - degreeV, vSpanIndex, controlPointsOfDerivative);

	for (int k = 0; k <= du; k++)
	{
		int dd = std::min(static_cast<int>(derivative - k), dv);
		for (int l = 0; l <= dd; l++)
		{
			derivatives[k][l] = XYZ();
			for (unsigned int i = 0; i <= degreeV - l; i++)
			{
				XYZ temp = XYZ();
				for (unsigned int j = 0; j <= degreeU - k; j++)
				{
					temp += allBasisFunctionsU[j][degreeU - k] * controlPointsOfDerivative[k][l][j][i];
				}
				derivatives[k][l] += allBasisFunctionsV[i][degreeV - l] * temp;
			}
		}
	}
}
