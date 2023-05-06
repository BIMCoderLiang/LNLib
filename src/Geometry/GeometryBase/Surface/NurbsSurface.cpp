#include "NurbsSurface.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "BsplineSurface.h"

void LNLib::NurbsSurface::GetPointOnSurface(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<std::vector<double>>& weights, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, XYZ& point)
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

	std::vector<std::vector<XYZW>> weightedPoints;
	for (unsigned int i = 0; i < controlPoints.size(); i++)
	{
		for (unsigned int j = 0; j < controlPoints[i].size(); j++)
		{
			weightedPoints[i][j] = XYZW(controlPoints[i][j], weights[i][j]);
		}
	}



	std::vector<XYZW> temp;
	for (unsigned int l = 0; l <= degreeV; l++)
	{
		temp[l] = XYZW(0, 0, 0, 0);
		for (unsigned int k = 0; k <= degreeU; k++)
		{
			temp[l] += basisFunctionsU[k] * weightedPoints[uSpanIndex - degreeU + k][vSpanIndex - degreeV + l];
		}
	}

	XYZW Sw = XYZW(0, 0, 0, 0);
	for (unsigned int l = 0; l <= degreeV; l++)
	{
		Sw += basisFunctionsV[l] * temp[l];
	}
	point = Sw.ToXYZ(true);
}


void LNLib::NurbsSurface::ComputeRationalSurfaceDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<std::vector<double>>& weights, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives)
{
	std::vector<std::vector<XYZW>> weightedPoints;
	weightedPoints.resize(controlPoints.size());
	for (unsigned int i = 0; i < controlPoints.size(); i++)
	{
		for (unsigned int j = 0; j < controlPoints[i].size(); j++)
		{
			weightedPoints[i][j] = XYZW(controlPoints[i][j], weights[i][j]);
		}
	}


	std::vector<std::vector<XYZW>> ders;
	ders.resize(degreeU + 1);
	for (unsigned int i = 0; i <= degreeU; i++)
	{
		ders[i].resize(degreeV + 1);
	}
	BsplineSurface::ComputeDerivatives(weightedPoints, knotVectorU, knotVectorV, degreeU, degreeV, paramU, paramV, derivative, ders);

	std::vector<std::vector<XYZ>> Aders;
	for (unsigned int i = 0; i <= degreeU; i++)
	{
		for (unsigned int j = 0; j <= degreeV; j++)
		{
			Aders[i][j] = ders[i][j].ToXYZ(false);
		}
	}

	std::vector<std::vector<double>> wders;
	for (unsigned int i = 0; i <= degreeU; i++)
	{
		for (unsigned int j = 0; j <= degreeV; j++)
		{
			wders[i][j] = ders[i][j].GetW();
		}
	}

	for (unsigned int k = 0; k <= derivative; k++)
	{
		for (unsigned int l = 0; l <= derivative - k; l++)
		{
			XYZ v = Aders[k][l];
			for (unsigned int j = 1; j <= l; j++)
			{
				v = v - MathUtils::Binomial(l, j) * wders[0][j] * derivatives[k][l - j];
			}

			for (unsigned int i = 1; i <= k; i++)
			{
				v = v - MathUtils::Binomial(k, i) * wders[i][0] * derivatives[k - i][l];

				XYZ v2 = XYZ();
				for (unsigned int j = 1; j <= l; j++)
				{
					v2 = v2 + MathUtils::Binomial(l, j) * wders[i][j] * derivatives[k - i][l - j];
				}

				v = v - MathUtils::Binomial(k, i) * v2;
			}

			derivatives[k][l] = v / wders[0][0];
		}
	}
}
