
#include "BezierCurve.h"
#include "ValidationUtils.h"
#include "XYZ.h"
#include "Polynomials.h"
#include "ValidationUtils.h"
#include <vector>

using namespace LNLib;

void BezierCurve::GetPointOnCurveByBernstein(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT, XYZ& point)
{
	std::vector<double> bernsteinArray;
	bernsteinArray.resize(controlPoints.size());
	Polynomials::AllBernstein(degree, paramT, bernsteinArray);

	XYZ temp(0,0,0);
	for (unsigned int k = 0; k <= degree; k++)
	{
		temp += bernsteinArray[k] * controlPoints[k];
	}
	point = temp;
}

void BezierCurve::GetPointOnCurveByDeCasteljau(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT, XYZ& point)
{
	std::vector<XYZ> temp;
	temp.resize(controlPoints.size());
	for (unsigned int i = 0; i <= degree; i++)
	{
		temp[i] = controlPoints[i];
	}
	for (unsigned int k = 1; k <= degree; k++)
	{
		for (unsigned int i = 0; i <= degree - k; i++)
		{
			temp[i] = (1.0 - paramT) * temp[i] + paramT * temp[i + 1];
		}
	}
	point = temp[0];
}
