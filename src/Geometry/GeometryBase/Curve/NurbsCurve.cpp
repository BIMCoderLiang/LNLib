#include "NurbsCurve.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "BsplineCurve.h"
#include <vector>

void LNLib::NurbsCurve::GetPointOnCurve(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double paramT, XYZ& point)
{
	int controlPointsSize = static_cast<int>(controlPoints.size());
	int n = controlPointsSize - 1;
	int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);
	std::vector<double> basisFunctions;
	basisFunctions.resize(degree + 1);
	Polynomials::BasisFunctions(spanIndex, degree, paramT, knotVector, basisFunctions);
	

	XYZW temp = XYZW(0,0,0,0);
	for (unsigned int j = 0; j <= degree; j++)
	{
		temp += basisFunctions[j] * controlPoints[j];
	}
	point = temp.ToXYZ(true);

}

void LNLib::NurbsCurve::ComputeRationalCurveDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives)
{
	int size = static_cast<int>(controlPoints.size());

	std::vector<XYZW> ders;
	ders.resize(degree + 1);
	BsplineCurve::ComputeDerivatives(degree, knotVector, controlPoints, paramT, derivative, ders);

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

void LNLib::NurbsCurve::InsertKnot(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double insertKnot, unsigned int times, std::vector<double>& insertedKnotVector, std::vector<XYZW>& updatedControlPoints)
{

	int np = static_cast<int>(controlPoints.size() - 1);
	int mp = np + degree + 1;

	int knotSpanIndex = Polynomials::GetKnotSpanIndex(np, degree, insertKnot, knotVector);
	int originMultiplicity = Polynomials::GetKnotMultiplicity(degree, insertKnot, knotVector);


	for (int i = 0; i <= knotSpanIndex; i++)
	{
		insertedKnotVector[i] = knotVector[i];
	}
	for (unsigned int i = 1; i <= times; i++)
	{
		insertedKnotVector[knotSpanIndex + i] = insertKnot;
	}
	for (int i = knotSpanIndex + 1; i <= mp; i++)
	{
		insertedKnotVector[i+times] = knotVector[i];
	}

	for (unsigned int i = 0; i <= knotSpanIndex - degree; i++)
	{
		updatedControlPoints[i] = controlPoints[i];
	}
	for (int i = knotSpanIndex - originMultiplicity; i <= np; i++)
	{
		updatedControlPoints[i+times] = controlPoints[i];
	}

	std::vector<XYZW> temp;
	temp.resize(degree - originMultiplicity + 1);
	for (unsigned int i = 0; i <= degree - originMultiplicity; i++)
	{
		temp[i] = controlPoints[knotSpanIndex - degree + i];
	}
	for (unsigned int j = 1; j <= times; j++)
	{
		int l = knotSpanIndex - degree + j;
		for (unsigned int i = 0; i <= degree - j - originMultiplicity; i++)
		{
			double alpha = (insertKnot - knotVector[l + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[l + i]);
			temp[i] = alpha * temp[i + 1] + (1 - alpha) * temp[i];
		}
		updatedControlPoints[l] = temp[0];
		updatedControlPoints[knotSpanIndex + times - j - originMultiplicity] = temp[degree - j - originMultiplicity];
	}

	int L = knotSpanIndex - degree + times;
	for (int i = L +1 ; i < knotSpanIndex - originMultiplicity; i++)
	{
		updatedControlPoints[i] = temp[i - L];
	}
}

void LNLib::NurbsCurve::GetPointOnCurveByInsertKnot(unsigned int degree,const std::vector<double>& knotVector, std::vector<XYZW>& controlPoints, double paramT, XYZ& point)
{
	int n = static_cast<int>(controlPoints.size() - 1);

	if (MathUtils::IsAlmostEqualTo(paramT, knotVector[0]))
	{
		point = controlPoints[0].ToXYZ(true);
		return;
	}
	if (MathUtils::IsAlmostEqualTo(paramT, knotVector[n + degree + 1]))
	{
		point = controlPoints[n].ToXYZ(true);
		return;
	}

	int knotSpanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);
	int originMultiplicity = Polynomials::GetKnotMultiplicity(degree, paramT, knotVector);

	int r = degree - originMultiplicity;
	std::vector<XYZW> temp;
	for (int i = 0; i <= r; i++)
	{
		temp[i] = controlPoints[knotSpanIndex - degree + i];
	}
	for (int j = 1; j <= r; j++)
	{
		for (int i = 0; i < r - j; i++)
		{
			double alpha = (paramT - knotVector[knotSpanIndex - degree + j + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[knotSpanIndex - degree + j + i]);
			temp[i] = alpha * temp[i + 1] + (1 - alpha) * temp[i];
		}
	}
	point = temp[0].ToXYZ(true);
}
