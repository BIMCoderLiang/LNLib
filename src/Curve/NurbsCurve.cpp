#include "NurbsCurve.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "BsplineCurve.h"
#include <vector>

void LNLib::NurbsCurve::GetPointOnCurve(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, const std::vector<double> weights, double paramT, XYZ& point)
{
	int controlPointsSize = static_cast<int>(controlPoints.size());
	int n = controlPointsSize - 1;
	int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);
	std::vector<double> basisFunctions;
	basisFunctions.resize(degree + 1);
	Polynomials::BasisFunctions(spanIndex, degree, paramT, knotVector, basisFunctions);

	std::vector<XYZW> weightedPoints;
	weightedPoints.resize(controlPointsSize);
	for (int i = 0; i < controlPointsSize; i++)
	{
		double w = weights[i];
		XYZ wc = controlPoints[i] * w;
		weightedPoints[i] = XYZW(wc, w);
	}

	XYZW temp = XYZW(0,0,0,0);
	for (unsigned int j = 0; j <= degree; j++)
	{
		temp += basisFunctions[j] * weightedPoints[j];
	}
	point = temp.ToXYZ(true);

}

void LNLib::NurbsCurve::ComputeRationalCurveDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, const std::vector<double> weights, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives)
{
	int size = static_cast<int>(controlPoints.size());

	std::vector<XYZW> weightedPoints;
	weightedPoints.resize(size);
	for (int i = 0; i < size; i++)
	{
		double w = weights[i];
		XYZ wc = controlPoints[i] * w;
		weightedPoints[i] = XYZW(wc, w);
	}

	std::vector<XYZW> ders;
	ders.resize(degree + 1);
	BsplineCurve::ComputeDerivatives(degree, knotVector, weightedPoints, paramT, derivative, ders);

	std::vector<XYZ> Aders;
	Aders.resize(degree + 1);
	for (unsigned int i = 0; i < ders.size(); i++)
	{
		Aders[i] = ders[i].ToXYZ(false);
	}
	std::vector<double> wders;
	wders.resize(degree + 1);
	for (unsigned int i = 0; i < ders.size(); i++)
	{
		wders[i] = ders[i].GetW();
	}

	for (unsigned int k = 0; k <= derivative; k++)
	{
		XYZ v = Aders[k];
		for (unsigned int i = 1; i <= k; i++)
		{
			v -=  MathUtils::Binomial(k, i) * wders[i] * derivatives[k - i];
		}
		derivatives[k] = v/wders[0];
	}

}
