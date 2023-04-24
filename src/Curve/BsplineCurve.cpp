#include "BsplineCurve.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "ValidationUtils.h"
#include <math.h>
#include <algorithm>

void LNLib::BsplineCurve::GetPointOnCurve(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT, const std::vector<double>& knotVector, XYZ& point)
{
	int n = static_cast<int>(controlPoints.size() - 1);
	int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);

	std::vector<double> basisFunctions;
	basisFunctions.resize(degree + 1);
	Polynomials::BasisFunctions(spanIndex, degree, paramT, knotVector, basisFunctions);

	for (unsigned int i = 0; i <= degree; i++)
	{
		point += basisFunctions[i] * controlPoints[spanIndex - degree + i];
	}
}

void LNLib::BsplineCurve::ComputeDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives)
{
	int du = std::min(derivative, degree);
	for (unsigned int k = degree + 1; k <= derivative; k++)
	{
		derivatives[k] = XYZ(0,0,0);
	}

	int n = static_cast<int>(controlPoints.size() - 1);
	int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);

	std::vector<std::vector<double>> derivativesOfBasisFunctions;
	derivativesOfBasisFunctions.resize(derivative + 1);
	for (unsigned i = 0; i <= derivative; i++)
	{
		derivativesOfBasisFunctions[i].resize(degree + 1);
	}
	Polynomials::BasisFunctionsDerivatives(spanIndex, degree, paramT,  du, knotVector, derivativesOfBasisFunctions);

	for (int k = 0; k <= du; k++)
	{
		derivatives[k] = XYZ(0, 0, 0);
		for (unsigned int j = 0; j <= degree; j++)
		{
			derivatives[k] += derivativesOfBasisFunctions[k][j] * controlPoints[spanIndex - degree + j];
		}
	}
}

void LNLib::BsplineCurve::ComputeDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZW>& derivatives)
{
	int du = std::min(derivative, degree);
	for (unsigned int k = degree + 1; k <= derivative; k++)
	{
		derivatives[k] = XYZW(0, 0, 0, 0);
	}

	int n = static_cast<int>(controlPoints.size() - 1);
	int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);

	std::vector<std::vector<double>> derivativesOfBasisFunctions;
	derivativesOfBasisFunctions.resize(derivative + 1);
	for (unsigned i = 0; i <= derivative; i++)
	{
		derivativesOfBasisFunctions[i].resize(degree + 1);
	}
	Polynomials::BasisFunctionsDerivatives(spanIndex, degree, paramT, du, knotVector, derivativesOfBasisFunctions);

	for (int k = 0; k <= du; k++)
	{
		derivatives[k] = XYZW(0, 0, 0, 0);
		for (unsigned int j = 0; j <= degree; j++)
		{
			derivatives[k] += derivativesOfBasisFunctions[k][j] * controlPoints[spanIndex - degree + j];
		}
	}
}

void LNLib::BsplineCurve::ComputeControlPointsOfDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, unsigned int derivative, unsigned int derMin, unsigned int derMax, std::vector<std::vector<XYZ>>& controlPointsOfDerivative)
{
	int range = derMax - derMin;
	for (int i = 0; i <= range; i++)
	{
		controlPointsOfDerivative[0][i] = controlPoints[derMin + i];
	}
	for (unsigned int k = 1; k <= derivative; k++)
	{
		int temp = degree - k + 1;
		for (unsigned int i = 0; i <= range - k; i++)
		{
			controlPointsOfDerivative[k][i] = temp * (controlPointsOfDerivative[k - 1][i + 1] - controlPointsOfDerivative[k - 1][i]) / (knotVector[derMin + i + degree + 1] - knotVector[derMin + i + k]);
		}
	}
}

void LNLib::BsplineCurve::ComputeDerivativesByAllBasisFunctions(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives)
{
	int du = std::min(derivative, degree);
	for (unsigned int k = degree + 1; k <= derivative; k++)
	{
		derivatives[k] = XYZ(0, 0, 0);
	}
	int n = static_cast<int>(controlPoints.size() - 1);
	int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);

	std::vector<std::vector<double>> allBasisFunctions;
	for (unsigned int i = 0; i <= degree; i++)
	{
		std::vector<double> basisFunctions;
		Polynomials::BasisFunctions(spanIndex, i, paramT, knotVector, basisFunctions);
		allBasisFunctions[i] = basisFunctions;
	}
	
	std::vector<std::vector<XYZ>> controlPointsOfDerivative;
	ComputeControlPointsOfDerivatives(degree, knotVector, controlPoints, du, spanIndex - degree, spanIndex, controlPointsOfDerivative);

	for (int k = 0; k <= du; k++)
	{
		derivatives[k] = XYZ(0, 0, 0);
		for (unsigned int j = 0; j <= degree - k; j++)
		{
			derivatives[k] += allBasisFunctions[j][degree - k] * controlPointsOfDerivative[k][j];
		}
	}
}
