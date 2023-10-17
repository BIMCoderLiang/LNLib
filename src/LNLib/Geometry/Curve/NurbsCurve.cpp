/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "NurbsCurve.h"
#include "Constants.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "Matrix4d.h"
#include "MathUtils.h"
#include "BezierCurve.h"
#include "BsplineCurve.h"
#include "Intersection.h"
#include "Projection.h"
#include "ValidationUtils.h"
#include "Interpolation.h"
#include "LNLibExceptions.h"
#include <vector>
#include <algorithm>

namespace LNLib
{
	double GetNode(int degree, const std::vector<double>& knotVector, int lastIndex)
	{
		double t = 0.0;
		for (int i = 0; i <= lastIndex; i++)
		{
			double sum = 0.0;
			for (int j = 1; j <= degree; j++)
			{
				sum += knotVector[i + j];
			}
			t = sum * (1.0 / degree);
		}
		return t;
	}
}

LNLib::XYZ LNLib::NurbsCurve::GetPointOnCurve(int degree, const std::vector<double>& knotVector, double paramT, const std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	XYZW weightPoint = BsplineCurve::GetPointOnCurve(degree, knotVector, paramT, controlPoints);
	return weightPoint.ToXYZ(true);
}

std::vector<LNLib::XYZ> LNLib::NurbsCurve::ComputeRationalCurveDerivatives(int degree, int derivative, const std::vector<double>& knotVector, double paramT, const std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(derivative > 0, "derivative", "derivative must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<LNLib::XYZ> derivatives(derivative + 1);
	std::vector<XYZW> ders = BsplineCurve::ComputeDerivatives(degree, derivative, knotVector, paramT, controlPoints);

	std::vector<XYZ> Aders(derivative + 1);
	for (int i = 0; i < ders.size(); i++)
	{
		Aders[i] = ders[i].ToXYZ(false);
	}
	std::vector<double> wders(derivative + 1);
	for (int i = 0; i < ders.size(); i++)
	{
		wders[i] = ders[i].GetW();
	}

	for (int k = 0; k <= derivative; k++)
	{
		XYZ v = Aders[k];
		for (int i = 1; i <= k; i++)
		{
			v = v - MathUtils::Binomial(k, i) * wders[i] * derivatives[k - i];
		}
		derivatives[k] = v / wders[0];
	}
	return derivatives;
}

void LNLib::NurbsCurve::InsertKnot(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double insertKnot, int times, std::vector<double>& insertedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(times > 0, "times", "Times must greater than zero.");

	int knotSpanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnot);
	int originMultiplicity = Polynomials::GetKnotMultiplicity(knotVector, insertKnot);

	if (originMultiplicity + times > degree)
	{
		times = degree - originMultiplicity;
	}

	insertedKnotVector.resize(knotVector.size() + times);
	for (int i = 0; i <= knotSpanIndex; i++)
	{
		insertedKnotVector[i] = knotVector[i];
	}
	for (int i = 1; i <= times; i++)
	{
		insertedKnotVector[knotSpanIndex + i] = insertKnot;
	}
	for (int i = knotSpanIndex + 1; i < knotVector.size(); i++)
	{
		insertedKnotVector[i + times] = knotVector[i];
	}

	updatedControlPoints.resize(controlPoints.size() + times);
	for (int i = 0; i <= knotSpanIndex - degree; i++)
	{
		updatedControlPoints[i] = controlPoints[i];
	}
	for (int i = knotSpanIndex - originMultiplicity; i < controlPoints.size(); i++)
	{
		updatedControlPoints[i + times] = controlPoints[i];
	}

	std::vector<XYZW> temp(degree - originMultiplicity + 1);
	for (int i = 0; i <= degree - originMultiplicity; i++)
	{
		temp[i] = controlPoints[knotSpanIndex - degree + i];
	}

	int L = 0;
	for (int j = 1; j <= times; j++)
	{
		L = knotSpanIndex - degree + j;
		for (int i = 0; i <= degree - j - originMultiplicity; i++)
		{
			double alpha = (insertKnot - knotVector[L + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
			temp[i] = alpha * temp[i + 1] + (1.0 - alpha) * temp[i];
		}
		updatedControlPoints[L] = temp[0];
		updatedControlPoints[knotSpanIndex + times - j - originMultiplicity] = temp[degree - j - originMultiplicity];
	}

	for (int i = L + 1; i < knotSpanIndex - originMultiplicity; i++)
	{
		updatedControlPoints[i] = temp[i - L];
	}
}

LNLib::XYZ LNLib::NurbsCurve::GetPointOnCurveByCornerCut(int degree, const std::vector<double>& knotVector, double paramT, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	if (MathUtils::IsAlmostEqualTo(paramT, knotVector[0]))
	{
		return controlPoints[0].ToXYZ(true);
	}
	int n = static_cast<int>(controlPoints.size() - 1);
	if (MathUtils::IsAlmostEqualTo(paramT, knotVector[n + degree + 1]))
	{
		return controlPoints[n].ToXYZ(true);
	}

	int knotSpanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, paramT);
	int originMultiplicity = Polynomials::GetKnotMultiplicity(knotVector, paramT);

	int times = degree - originMultiplicity;
	std::vector<XYZW> temp(times + 1);
	for (int i = 0; i <= times; i++)
	{
		temp[i] = controlPoints[knotSpanIndex - degree + i];
	}
	for (int j = 1; j <= times; j++)
	{
		for (int i = 0; i <= times - j; i++)
		{
			double alpha = (paramT - knotVector[knotSpanIndex - degree + j + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[knotSpanIndex - degree + j + i]);
			temp[i] = alpha * temp[i + 1] + (1.0 - alpha) * temp[i];
		}
	}
	return temp[0].ToXYZ(true);
}

void LNLib::NurbsCurve::RefineKnotVector(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& insertKnotElements, std::vector<double>& insertedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(insertKnotElements.size() > 0, "insertKnotElements", "insertKnotElements size must greater than zero.");

	int n = controlPoints.size() - 1;
	int m = n + degree + 1;
	int r = insertKnotElements.size() - 1;

	int a = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnotElements[0]);
	int b = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnotElements[r]) + 1;

	insertedKnotVector.resize(m + r + 2);
	for (int j = 0; j <= a; j++)
	{
		insertedKnotVector[j] = knotVector[j];
	}
	for (int j = b + degree; j <= m; j++)
	{
		insertedKnotVector[j + r + 1] = knotVector[j];
	}

	updatedControlPoints.resize(n + r + 2);
	for (int j = 0; j <= a - degree; j++)
	{
		updatedControlPoints[j] = controlPoints[j];
	}
	for (int j = b - 1; j <= n; j++)
	{
		updatedControlPoints[j + r + 1] = controlPoints[j];
	}

	int i = b + degree - 1;
	int k = b + degree + r;
	for (int j = r; j >= 0; j--)
	{
		while (insertKnotElements[j] <= knotVector[i] && i > a)
		{
			updatedControlPoints[k - degree - 1] = controlPoints[i - degree - 1];
			insertedKnotVector[k] = knotVector[i];
			k = k - 1;
			i = i - 1;
		}

		updatedControlPoints[k - degree - 1] = updatedControlPoints[k - degree];
		for (int l = 1; l <= degree; l++)
		{
			int ind = k - degree + l;
			double alpha = insertedKnotVector[k + l] - insertKnotElements[j];
			if (MathUtils::IsAlmostEqualTo(abs(alpha), 0.0))
			{
				updatedControlPoints[ind - 1] = updatedControlPoints[ind];
			}
			else
			{
				alpha = alpha / (insertedKnotVector[k + l] - knotVector[i - degree + l]);
				updatedControlPoints[ind - 1] = alpha * updatedControlPoints[ind - 1] + (1.0 - alpha) * updatedControlPoints[ind];
			}
		}
		insertedKnotVector[k] = insertKnotElements[j];
		k = k - 1;
	}
}

std::vector<std::vector<LNLib::XYZW>> LNLib::NurbsCurve::DecomposeToBeziers(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<std::vector<XYZW>> decomposedControlPoints(controlPoints.size() - degree, std::vector<XYZW>(degree + 1));

	int n = controlPoints.size() - 1;
	int m = n + degree + 1;

	int a = degree;
	int b = degree + 1;

	int nb = 0;
	for (int i = 0; i <= degree; i++)
	{
		decomposedControlPoints[nb][i] = controlPoints[i];
	}

	while (b < m)
	{
		int i = b;
		while (b < m && MathUtils::IsAlmostEqualTo(knotVector[b + 1], knotVector[b]))
		{
			b++;
		}
		int multi = b - i + 1;
		if (multi < degree)
		{
			double numerator = knotVector[b] - knotVector[a];
			std::vector<double> alphaVector(degree + 1);
			for (int j = degree; j > multi; j--)
			{
				alphaVector[j - multi - 1] = numerator / (knotVector[a + j] - knotVector[a]);
			}

			int r = degree - multi;
			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				int s = multi + j;
				for (int k = degree; k >= s; k--)
				{
					double alpha = alphaVector[k - s];
					decomposedControlPoints[nb][k] = alpha * decomposedControlPoints[nb][k] + (1.0 - alpha) * decomposedControlPoints[nb][k - 1];
				}

				if (b < m)
				{
					decomposedControlPoints[nb + 1][save] = decomposedControlPoints[nb][degree];
				}
			}

			nb++;
			if (b < m)
			{
				for (int i = degree - multi; i <= degree; i++)
				{
					decomposedControlPoints[nb][i] = controlPoints[b - degree + i];
				}

				a = b;
				b += 1;
			}
		}
	}
	return decomposedControlPoints;
}

void LNLib::NurbsCurve::RemoveKnot(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double removeKnot, int times, std::vector<double>& restKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT_RANGE(removeKnot, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(times > 0, "times", "Times must greater than zero.");

	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);
	int n = static_cast<int>(controlPoints.size() - 1);

	int order = degree + 1;
	int s = Polynomials::GetKnotMultiplicity(knotVector, removeKnot);
	int r = Polynomials::GetKnotSpanIndex(degree, knotVector, removeKnot);

	int first = r - degree;
	int last = r - s;

	restKnotVector = knotVector;
	int m = n + degree + 1;
	for (int k = r + 1; k <= m; k++)
	{
		restKnotVector[k - times] = restKnotVector[k];
	}
	for (int i = 0; i < times; i++)
	{
		restKnotVector.pop_back();
	}

	updatedControlPoints = controlPoints;
	std::vector<XYZW> temp(2 * degree + 1);

	int t = 0;
	for (t = 0; t < times; t++)
	{
		int off = first - 1;
		temp[0] = controlPoints[off];
		temp[last + 1 - off] = controlPoints[last + 1];
		int i = first;
		int j = last;
		int ii = 1;
		int jj = last - off;
		bool remflag = false;

		while (j - i >= t)
		{
			double alphai = (removeKnot - knotVector[i]) / (knotVector[i + order + t] - knotVector[i]);
			double alphaj = (removeKnot - knotVector[j - t]) / (knotVector[j + order] - knotVector[j - t]);

			temp[ii] = (controlPoints[i] - (1.0 - alphai) * temp[ii - 1]) / alphai;
			temp[jj] = (controlPoints[j] - alphaj * temp[jj + 1]) / (1.0 - alphaj);

			i = i + 1;
			ii = ii + 1;

			j = j - 1;
			jj = jj - 1;
		}

		if (j - i < t)
		{
			if (MathUtils::IsLessThanOrEqual(temp[ii - 1].Distance(temp[jj + 1]), tol))
			{
				remflag = true;
			}
		}
		else
		{
			double alphai = (removeKnot - knotVector[i]) / (knotVector[i + order + t] - knotVector[i]);
			if (MathUtils::IsLessThanOrEqual(controlPoints[i].Distance(alphai * temp[ii + t + 1] + (1.0 * alphai) * temp[ii - 1]), tol))
			{
				remflag = true;
			}
		}

		if (remflag)
		{
			i = first;
			j = last;

			while (j - i > t)
			{
				updatedControlPoints[i] = temp[i - off];
				updatedControlPoints[j] = temp[j - off];
				i = i + 1;
				j = j - 1;
			}
		}
		first = first - 1;
		last = last + 1;
	}

	if (t == 0)
	{
		return;
	}

	int j = (2 * r - s - degree) / 2;
	int i = j;

	for (int k = 1; k < t; k++)
	{
		if (k % 2 == 1)
		{
			i = i + 1;
		}
		else
		{
			j = j - 1;
		}
	}

	for (int k = i + 1; k <= n; k++)
	{
		updatedControlPoints[j] = controlPoints[k];
		j = j + 1;
	}
	for (int i = 0; i < t; i++)
	{
		updatedControlPoints.pop_back();
	}
	return;
}

void LNLib::NurbsCurve::ElevateDegree(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, int times, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(times > 0, "times", "Times must greater than zero.");

	int n = controlPoints.size() - 1;
	int m = n + degree + 1;
	int ph = degree + times;
	int ph2 = floor(ph / 2);

	std::vector<std::vector<double>> bezalfs(degree + times + 1, std::vector<double>(degree + 1));
	bezalfs[0][0] = bezalfs[ph][degree] = 1.0;

	for (int i = 1; i <= ph2; i++)
	{
		double inv = 1.0 / MathUtils::Binomial(ph, i);
		int mpi = std::min(degree, i);

		for (int j = std::max(0, i - times); j <= mpi; j++)
		{
			bezalfs[i][j] = inv * MathUtils::Binomial(degree, j) * MathUtils::Binomial(times, i - j);
		}
	}

	for (int i = ph2 + 1; i <= ph - 1; i++)
	{
		int mpi = std::min(degree, i);
		for (int j = std::max(0, i - times); j <= mpi; j++)
		{
			bezalfs[i][j] = bezalfs[ph - i][degree - j];
		}
	}

	int mh = ph;
	int kind = ph + 1;
	int r = -1;
	int a = degree;
	int b = degree + 1;
	int cind = 1;
	double ua = knotVector[0];

	updatedControlPoints.resize(n + n * times);
	updatedControlPoints[0] = controlPoints[0];

	updatedKnotVector.resize(n + n * times + ph + 1);
	for (int i = 0; i <= ph; i++)
	{
		updatedKnotVector[i] = ua;
	}

	std::vector<XYZW> bpts(degree + 1);
	for (int i = 0; i <= degree; i++)
	{
		bpts[i] = controlPoints[i];
	}

	std::vector<XYZW> nextbpts(degree - 1);

	while (b < m)
	{
		int i = b;
		while (b < m && MathUtils::IsAlmostEqualTo(knotVector[b], knotVector[b + 1]))
		{
			b = b + 1;
		}
		int mul = b - i + 1;
		mh += mul + times;
		double ub = knotVector[b];

		int oldr = r;
		r = degree - mul;

		int lbz = oldr > 0 ? floor((oldr + 2) / 2) : 1;
		int rbz = r > 0 ? floor(ph - (r + 1) / 2) : ph;

		if (r > 0)
		{
			double numer = ub - ua;
			std::vector<double> alfs(degree - 1);
			for (int k = degree; k > mul; k--)
			{
				alfs[k - mul - 1] = numer / (knotVector[a + k] - ua);
			}
			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				int s = mul + j;

				for (int k = degree; k >= s; k--)
				{
					bpts[k] = alfs[k - s] * bpts[k] + (1.0 - alfs[k - s]) * bpts[k - 1];
				}
				nextbpts[save] = bpts[degree];
			}
		}

		std::vector<XYZW> ebpts(degree + times + 1);
		for (int i = lbz; i <= ph; i++)
		{
			ebpts[i] = XYZW(0.0, 0.0, 0.0, 0.0);
			int mpi = std::min(degree, i);
			for (int j = std::max(0, i - times); j <= mpi; j++)
			{
				ebpts[i] += bezalfs[i][j] * bpts[j];
			}
		}

		if (oldr > 1)
		{
			int first = kind - 2;
			int last = kind;
			double den = ub - ua;
			double bet = (ub - updatedKnotVector[kind - 1]) / den;

			for (int tr = 1; tr < oldr; tr++)
			{
				int i = first;
				int j = last;
				int kj = j - kind + 1;

				while (j - i > tr)
				{
					if (i < cind)
					{
						double alf = (ub - updatedKnotVector[i]) / (ua - updatedKnotVector[i]);
						updatedControlPoints[i] = alf * updatedControlPoints[i] + (1.0 - alf) * updatedControlPoints[i - 1];
					}

					if (j >= lbz)
					{
						if (j - tr <= kind - ph + oldr)
						{
							double gam = (ub - updatedKnotVector[j - tr]) / den;
							ebpts[kj] = gam * ebpts[kj] + (1.0 - gam) * ebpts[kj + 1];
						}
						else
						{
							ebpts[kj] = bet * ebpts[kj] + (1.0 - bet) * ebpts[kj + 1];
						}
					}

					i = i + 1;
					j = j - 1;
					kj = kj - 1;
				}

				first -= 1;
				last += 1;
			}
		}

		if (a != degree)
		{
			for (int i = 0; i < ph - oldr; i++)
			{
				updatedKnotVector[kind++] = ua;
			}
		}

		for (int j = lbz; j <= rbz; j++)
		{
			updatedControlPoints[cind++] = ebpts[j];
		}

		if (b < m)
		{
			for (int j = 0; j < r; j++)
			{
				bpts[j] = nextbpts[j];
			}
			for (int j = r; j <= degree; j++)
			{
				bpts[j] = controlPoints[b - degree + j];
			}

			a = b;
			b = b + 1;
			ua = ub;
		}
		else
		{
			for (int i = 0; i <= ph; i++)
			{
				updatedKnotVector[kind + i] = ub;
			}
		}
	}

	int diffcp = n + n * times - (mh - ph);
	int diffkv = n + n * times + ph + 1 - (mh + 1);
	for (int i = 0; i < diffcp; i++)
	{
		updatedControlPoints.pop_back();
	}
	for (int i = 0; i < diffkv; i++)
	{
		updatedKnotVector.pop_back();
	}
}

bool LNLib::NurbsCurve::ReduceDegree(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(ValidationUtils::IsValidDegreeReduction(degree), "degree", "Degree must greater than one.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);

	int size = controlPoints.size();
	bool isBezier = ValidationUtils::IsValidBezier(degree, size);
	if (!isBezier) return false;

	int r = floor((degree - 1) / 2);
	updatedControlPoints.resize(degree);
	updatedControlPoints[0] = controlPoints[0];
	updatedControlPoints[degree - 1] = controlPoints[degree];

	std::vector<double> alpha(degree);
	for (int i = 0; i < degree; i++)
	{
		alpha[i] = double(i) / double(degree);
	}
	double error = 0.0;
	if (degree % 2 == 0)
	{
		for (int i = 1; i < r + 1; i++)
		{
			updatedControlPoints[i] = (controlPoints[i] - alpha[i] * updatedControlPoints[i - 1]) / (1 - alpha[i]);
		}
		for (int i = degree - 2; i > r; i--)
		{
			updatedControlPoints[i] = (controlPoints[i + 1] - (1 - alpha[i+1]) * updatedControlPoints[i + 1]) / alpha[i + 1];
		}
		error = (controlPoints[r + 1].Distance(0.5 * (updatedControlPoints[r] + updatedControlPoints[r + 1])));
		double c = MathUtils::Binomial(degree, r + 1);
		error = error * (c * pow(0.5, r + 1) * pow(1 - 0.5, degree - r - 1));
	}
	else
	{
		for (int i = 1; i < r; i++)
		{
			updatedControlPoints[i] = (controlPoints[i] - alpha[i] * updatedControlPoints[i - 1]) / (1 - alpha[i]);
		}
		for (int i = degree - 2; i > r; i--)
		{
			updatedControlPoints[i] = (controlPoints[i + 1] - (1 - alpha[i + 1]) * updatedControlPoints[i + 1]) / alpha[i + 1];
		}
		XYZW PLr = (controlPoints[r] - alpha[r] * updatedControlPoints[r - 1]) / (1 - alpha[r]);
		XYZW PRr = (controlPoints[r + 1] - (1 - alpha[r + 1]) * updatedControlPoints[r + 1]) / alpha[r + 1];
		updatedControlPoints[r] = 0.5 * (PLr + PRr);
		error = PLr.Distance(PRr);
		double maxU = (degree - sqrt(degree)) / (2 * degree);
		error = error * 0.5 * (1 - alpha[r]) * (MathUtils::Binomial(degree, r) * pow(maxU, r) * pow(1 - maxU, r + 1) * (1 - 2 * maxU));
	}
	if (error > tol) return false;

	auto map = Polynomials::GetKnotMultiplicityMap(knotVector);
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		double u = it->first;
		int count = it->second - 1;
		for (int i = 0; i < count; i++)
		{
			updatedKnotVector.emplace_back(u);
		}
	}
	return true;
}

void LNLib::NurbsCurve::EquallyTessellate(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<XYZ>& tessellatedPoints, std::vector<double>& correspondingKnots)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<double> uniqueKv = knotVector;
	uniqueKv.erase(unique(uniqueKv.begin(), uniqueKv.end()), uniqueKv.end());
	int size = uniqueKv.size();
	int intervals = 100;
	for (int i = 0; i < size - 1; i++)
	{
		double currentU = uniqueKv[i];
		double nextU = uniqueKv[i + 1];
		double step = (nextU - currentU) / intervals;
		for (int j = 0; j < intervals; j++)
		{
			double u = currentU + step * j;
			correspondingKnots.emplace_back(u);
			tessellatedPoints.emplace_back(GetPointOnCurve(degree, knotVector, u, controlPoints));
		}
	}
	correspondingKnots.emplace_back(knotVector[knotVector.size() - 1]);
	tessellatedPoints.emplace_back(const_cast<XYZW&>(controlPoints[controlPoints.size() - 1]).ToXYZ(true));
}

double LNLib::NurbsCurve::GetParamOnCurve(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, const XYZ& givenPoint)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	double minValue = Constants::MaxDistance;

	int maxIterations = 10;
	double paramT = Constants::DoubleEpsilon;
	double minParam = knotVector[0];
	double maxParam = knotVector[knotVector.size() - 1];

	std::vector<XYZ> tessellatedPoints;
	std::vector<double> correspondingKnots;
	EquallyTessellate(degree, knotVector, controlPoints, tessellatedPoints, correspondingKnots);
	for (int i = 0; i < tessellatedPoints.size() - 1; i++)
	{
		double currentU = correspondingKnots[i];
		double nextU = correspondingKnots[i + 1];

		XYZ currentPoint = tessellatedPoints[i];
		XYZ nextPoint = tessellatedPoints[i + 1];

		XYZ vector1 = currentPoint - givenPoint;
		XYZ vector2 = nextPoint - currentPoint;
		double dot = vector1.DotProduct(vector2);

		XYZ projectPoint;
		double projectU;

		if (dot < 0)
		{
			projectPoint = currentPoint;
			projectU = currentU;
		}
		else if (dot > 1)
		{
			projectPoint = nextPoint;
			projectU = nextU;
		}
		else
		{
			projectPoint = currentPoint + dot * vector1.Normalize();
			projectU = currentU + (nextU - currentU) * dot;
		}

		double distance = (givenPoint - projectPoint).Length();
		if (distance < minValue)
		{
			minValue = distance;
			paramT = projectU;
		}
	}

	bool isClosed = ValidationUtils::IsClosed(controlPoints);
	double a = minParam;
	double b = maxParam;

	int counters = 0;
	while (counters < maxIterations)
	{
		std::vector<XYZ> derivatives = ComputeRationalCurveDerivatives(degree, 2, knotVector, paramT, controlPoints);
		XYZ difference = derivatives[0] - givenPoint;
		double f = derivatives[1].DotProduct(difference);

		double condition1 = difference.Length();
		double condition2 = std::abs(f / (derivatives[1].Length() * condition1));

		if (condition1 < Constants::DistanceEpsilon &&
			condition2 < Constants::DistanceEpsilon)
		{
			return paramT;
		}

		double df = derivatives[2].DotProduct(difference) + derivatives[1] * derivatives[1];
		double temp = paramT - f / df;

		if (isClosed)
		{
			if (temp < a)
			{
				temp = a;
			}
			if (temp > b)
			{
				temp = b;
			}
		}
		else
		{
			if (temp < a)
			{
				temp = b - (a - temp);
			}
			if (temp > b)
			{
				temp = a + (temp - b);
			}
		}

		double condition4 = ((temp - paramT) * derivatives[1]).Length();
		if (condition4 < Constants::DistanceEpsilon) {
			return paramT;
		}

		paramT = temp;
		counters++;
	}
	return paramT;
}

void LNLib::NurbsCurve::CreateTransformed(const std::vector<XYZW>& controlPoints, const Matrix4d& matrix, std::vector<XYZW>& transformedControlPoints)
{
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");

	int size = static_cast<int>(controlPoints.size());
	transformedControlPoints.resize(size);

	Matrix4d tempMatrix = matrix;
	for (int i = 0; i < size; i++)
	{
		XYZW temp = controlPoints[i];
		transformedControlPoints[i] = XYZW(tempMatrix.OfPoint(temp.ToXYZ(true)), temp.GetW());
	}
}

void LNLib::NurbsCurve::Reparameterization(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double alpha, double beta, double gamma, double delta, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(alpha * delta, gamma * beta), "coefficient", "(alpha * delta - gamma * beta) must greater than zero");

	updatedKnotVector.resize(knotVector.size());
	for (int i = 0; i < knotVector.size(); i++)
	{
		updatedKnotVector[i] = (alpha * knotVector[i] + beta) / (gamma * knotVector[i] + delta);
	}

	updatedControlPoints.resize(controlPoints.size());
	for (int i = 0; i < controlPoints.size(); i++)
	{
		double temp = 1.0;
		for (int j = 1; j <= degree; j++)
		{
			double lambda = updatedKnotVector[i + j] * gamma - alpha;
			temp = temp * lambda;
		}
		double newW = abs(controlPoints[i].GetW() * temp);
		updatedControlPoints[i] = XYZW(const_cast<XYZW&>(controlPoints[i]).ToXYZ(true), newW);
	}
}

void LNLib::NurbsCurve::ReverseKnotVector(const std::vector<double>& knotVector, std::vector<double> reversedKnotVector)
{
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");

	int size = static_cast<int>(knotVector.size());

	double min = knotVector[0];
	double max = knotVector[size - 1];

	reversedKnotVector.resize(size);
	reversedKnotVector[0] = min;
	for (int i = 1; i < size; i++)
	{
		reversedKnotVector[i] = reversedKnotVector[i - 1] + (knotVector[size - i] - knotVector[size - i - 1]);
	}
}

void LNLib::NurbsCurve::ReverseControlPoints(const std::vector<XYZW>& controlPoints, std::vector<XYZW> reversedControlPoints)
{
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");

	reversedControlPoints = controlPoints;
	std::reverse(reversedControlPoints.begin(), reversedControlPoints.end());
}

void LNLib::NurbsCurve::Reverse(const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& reversedKnotVector, std::vector<XYZW>& reversedControlPoints)
{
	ReverseKnotVector(knotVector, reversedKnotVector);
	ReverseControlPoints(controlPoints, reversedControlPoints);
}

bool LNLib::NurbsCurve::CreateArc(const XYZ& center, const XYZ& xAxis, const XYZ& yAxis, double startRad, double endRad, double xRadius, double yRadius, int& degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(!xAxis.IsZero(), "xAxis", "xAxis must not be zero vector.");
	VALIDATE_ARGUMENT(!yAxis.IsZero(), "yAxis", "yAxis must not be zero vector.")
		VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(endRad, startRad), "endRad", "endRad must greater than startRad.");
	double theta = endRad - startRad;
	VALIDATE_ARGUMENT_RANGE(theta, 0, 2 * Constants::Pi);
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(xRadius, 0.0), "xRadius", "xRadius must greater than zero.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(yRadius, 0.0), "yRadius", "yRadius must greater than zero.");

	int narcs = 0;
	if (MathUtils::IsLessThanOrEqual(theta, Constants::Pi / 2.0))
	{
		narcs = 1;
	}
	else
	{
		if (MathUtils::IsLessThanOrEqual(theta, Constants::Pi))
		{
			narcs = 2;
		}
		else if (MathUtils::IsLessThanOrEqual(theta, 3 * Constants::Pi / 2.0))
		{
			narcs = 3;
		}
		else
		{
			narcs = 4;
		}
	}
	double dtheta = theta / (double)narcs;
	int n = 2 * narcs;

	degree = 2;
	knotVector.resize(n + degree + 1 + 1);
	controlPoints.resize(n + 1);

	double w1 = cos(dtheta / 2.0);
	XYZ nX = const_cast<XYZ&>(xAxis).Normalize();
	XYZ nY = const_cast<XYZ&>(yAxis).Normalize();
	XYZ P0 = center + xRadius * cos(startRad) * nX + yRadius * sin(startRad) * nY;
	XYZ T0 = -sin(startRad) * nX + cos(startRad) * nY;

	controlPoints[0] = XYZW(P0, 1);
	int index = 0;
	double angle = startRad;
	for (int i = 1; i <= narcs; i++)
	{
		angle += dtheta;
		XYZ P2 = center + xRadius * cos(angle) * nX + yRadius * sin(angle) * nY;
		controlPoints[index + 2] = XYZW(P2, 1);
		XYZ T2 = -sin(angle) * nX + cos(angle) * nY;

		double param0, param2 = 0.0;
		XYZ P1;
		CurveCurveIntersectionType type = Intersection::ComputeRays(P0, T0, P2, T2, param0, param2, P1);
		if (type != CurveCurveIntersectionType::Intersecting) return false;
		controlPoints[index + 1] = XYZW(P1, w1);
		index = index + 2;
		if (i < narcs)
		{
			P0 = P2;
			T0 = T2;
		}
	}

	int j = 2 * narcs + 1;

	for (int i = 0; i < 3; i++)
	{
		knotVector[i] = 0.0;
		knotVector[i + j] = 1.0;
	}

	switch (narcs)
	{
	case 1:
		break;
	case 2:
		knotVector[3] = knotVector[4] = 0.5;
		break;
	case 3:
		knotVector[3] = knotVector[4] = 1 / 3;
		knotVector[5] = knotVector[6] = 2 / 3;
		break;
	case 4:
		knotVector[3] = knotVector[4] = 0.25;
		knotVector[5] = knotVector[6] = 0.5;
		knotVector[7] = knotVector[8] = 0.75;
		break;
	}
	return true;
}

bool LNLib::NurbsCurve::CreateOneConicArc(const XYZ& start, const XYZ& startTangent, const XYZ& end, const XYZ& endTangent, const XYZ& pointOnConic, XYZ& projectPoint, double& projectPointWeight)
{
	VALIDATE_ARGUMENT(!startTangent.IsZero(), "startTangent", "StartTangent must not be zero vector.");
	VALIDATE_ARGUMENT(!endTangent.IsZero(), "endTangent", "EndTangent must not be zero vector.")

		double param0, param1 = 0.0;
	XYZ point = XYZ(0, 0, 0);
	CurveCurveIntersectionType type = Intersection::ComputeRays(start, startTangent, end, endTangent, param0, param1, point);

	XYZ pDiff = end - start;
	double alf0, alf2 = 0.0;
	XYZ dummy = XYZ(0, 0, 0);
	if (type == CurveCurveIntersectionType::Intersecting)
	{
		XYZ v1p = pointOnConic - point;
		type = Intersection::ComputeRays(point, v1p, start, pDiff, alf0, alf2, dummy);
		if (type == CurveCurveIntersectionType::Intersecting)
		{
			double a = sqrt(alf2 / (1.0 - alf2));
			double u = a / (1.0 + a);
			double num = (1.0 - u) * (1.0 - u) * (pointOnConic - start).DotProduct(point - pointOnConic) + u * u * (pointOnConic - end).DotProduct(point - pointOnConic);
			double den = 2.0 * u * (1.0 - u) * (point - pointOnConic).DotProduct(point - pointOnConic);
			projectPoint = point;
			projectPointWeight = num / den;
			return true;
		}
	}
	else if (type == CurveCurveIntersectionType::Parallel)
	{
		type = Intersection::ComputeRays(pointOnConic, startTangent, start, pDiff, alf0, alf2, dummy);
		if (type == CurveCurveIntersectionType::Intersecting)
		{
			double a = sqrt(alf2 / (1.0 - alf2));
			double u = a / (1.0 + a);
			double b = 2.0 * u * (1.0 - u);
			b = -alf0 * (1.0 - b) / b;
			projectPoint = b * startTangent;
			projectPointWeight = 0.0;
			return true;
		}
	}
	return false;
}

void LNLib::NurbsCurve::SplitArc(const XYZ& start, const XYZ& projectPoint, double projectPointWeight, const XYZ& end, XYZ& insertPointAtStartSide, XYZ& splitPoint, XYZ& insertPointAtEndSide, double insertWeight)
{
	insertPointAtStartSide = start + projectPoint;
	insertPointAtEndSide = end + projectPoint;
	splitPoint = (insertPointAtStartSide + insertPointAtEndSide) / 2.0;
	insertWeight = sqrt(1 + projectPointWeight) / 2.0;
}

bool LNLib::NurbsCurve::CreateOpenConic(const XYZ& start, const XYZ& startTangent, const XYZ& end, const XYZ& endTangent, const XYZ& pointOnConic, int& degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(!startTangent.IsZero(), "startTangent", "StartTangent must not be zero vector.");
	VALIDATE_ARGUMENT(!endTangent.IsZero(), "endTangent", "EndTangent must not be zero vector.")

		XYZ P1;
	double w1 = 0.0;
	bool isCreated = CreateOneConicArc(start, startTangent, end, endTangent, pointOnConic, P1, w1);
	if (!isCreated) return false;

	int nsegs = 0;
	if (MathUtils::IsLessThanOrEqual(w1, -1.0))
	{
		return false;
	}
	if (MathUtils::IsGreaterThanOrEqual(w1, 1.0))
	{
		nsegs = 1;
	}
	else
	{
		XYZ v1 = (P1 - start).Normalize();
		XYZ v2 = (end - P1).Normalize();
		double rad = v1.AngleTo(v2);
		if (MathUtils::IsGreaterThan(w1, 0.0) && rad > MathUtils::AngleToRadians(60))
		{
			nsegs = 1;
		}
		else if (MathUtils::IsLessThan(w1, 0.0) && rad > MathUtils::AngleToRadians(90))
		{
			nsegs = 4;
		}
		else
		{
			nsegs = 2;
		}
	}

	int n = 2 * nsegs;
	int j = 2 * nsegs + 1;

	degree = 2;
	knotVector.resize(j + degree + 1);
	controlPoints.resize(n + 1);

	for (int i = 0; i < 3; i++)
	{
		knotVector[i] = 0.0;
		knotVector[i + j] = 1.0;
	}

	controlPoints[0] = XYZW(start, 1.0);
	controlPoints[n] = XYZW(end, 1.0);

	if (nsegs == 1)
	{
		controlPoints[1] = XYZW(P1, w1);
		return true;
	}

	XYZ Q1, R1, S;
	double wqr = 0.0;
	SplitArc(start, P1, w1, end, Q1, S, R1, wqr);

	if (nsegs == 2)
	{
		controlPoints[2] = XYZW(S, 1.0);
		controlPoints[1] = XYZW(Q1, wqr);
		controlPoints[3] = XYZW(R1, wqr);

		knotVector[3] = knotVector[4] = 0.5;
		return true;
	}

	if (nsegs == 4)
	{
		controlPoints[4] = XYZW(S, 1.0);
		w1 = wqr;

		XYZ HQ1, HR1, HS;
		SplitArc(start, Q1, w1, S, HQ1, HS, HR1, wqr);
		controlPoints[2] = XYZW(HS, 1.0);
		controlPoints[1] = XYZW(HQ1, wqr);
		controlPoints[3] = XYZW(HR1, wqr);

		SplitArc(S, R1, w1, end, HQ1, HS, HR1, wqr);
		controlPoints[6] = XYZW(HS, 1.0);
		controlPoints[5] = XYZW(HQ1, wqr);
		controlPoints[7] = XYZW(HR1, wqr);

		for (int i = 0; i < 2; i++)
		{
			knotVector[i + 3] = 0.25;
			knotVector[i + 5] = 0.5;
			knotVector[i + 7] = 0.75;
		}
		return true;
	}
	return false;
}

void LNLib::NurbsCurve::GlobalInterpolation(int degree, const std::vector<XYZ>& throughPoints, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints.size() > degree, "throughPoints", "ThroughPoints size must greater than degree.");

	int size = throughPoints.size();
	int n = size - 1;

	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	knotVector = Interpolation::AverageKnotVector(degree, uk);

	std::vector<std::vector<double>> A(size, std::vector<double>(size));
	for (int i = 1; i < n; i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, uk[i]);
		std::vector<double> basis = Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i]);

		for (int j = 0; j <= degree; j++)
		{
			A[i][spanIndex - degree + j] = basis[j];
		}
	}
	A[0][0] = 1.0;
	A[n][n] = 1.0;

	std::vector<std::vector<double>> right(size, std::vector<double>(3));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			right[i][j] = throughPoints[i][j];
		}
	}

	controlPoints.resize(size);
	std::vector<std::vector<double>> result = MathUtils::SolveLinearSystem(A, right);
	for (int i = 0; i < result.size(); i++)
	{
		XYZ temp = XYZ(0, 0, 0);
		for (int j = 0; j < 3; j++)
		{
			temp[j] = result[i][j];
		}
		controlPoints[i] = XYZW(temp, 1.0);
	}
}

void LNLib::NurbsCurve::GlobalInterpolation(int degree, const std::vector<XYZ>& throughPoints, const std::vector<XYZ>& tangents, double tangentFactor, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(throughPoints.size() > degree, "throughPoints", "ThroughPoints size must greater than degree.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(tangentFactor, 0.0), "tangentFactor", "TangentFactor must greater than zero.");

	std::vector<XYZ> unitTangents(tangents.size());
	for (int i = 0; i < tangents.size(); i++)
	{
		unitTangents[i] = const_cast<XYZ&>(tangents[i]).Normalize();
	}

	int size = throughPoints.size();
	int n = 2 * size;

	controlPoints.resize(n);
	knotVector.resize(n + degree + 1);

	std::vector<double> uk(size);
	uk[size - 1] = 1.0;
	double d = Interpolation::GetTotalChordLength(throughPoints) * tangentFactor;
	for (int i = 1; i <= size - 2; i++)
	{
		uk[i] = uk[i - 1] + (throughPoints[i].Distance(throughPoints[i - 1])) / d;
	}
	switch (degree)
	{
	case 2:
	{
		for (int i = 0; i <= degree; i++)
		{
			knotVector[i] = 0.0;
			knotVector[knotVector.size() - 1 - i] = 1.0;
		}
		for (int i = 0; i < size - 1; i++)
		{
			knotVector[2 * i + degree] = uk[i];
			knotVector[2 * i + degree + 1] = (uk[i] + uk[i + 1]) / 2.0;
		}
		break;
	}

	case 3:
	{
		for (int i = 0; i <= degree; i++)
		{
			knotVector[i] = 0.0;
			knotVector[knotVector.size() - 1 - i] = 1.0;
		}
		for (int i = 1; i < size - 1; i++)
		{
			knotVector[degree + 2 * i] = (2 * uk[i] + uk[i + 1]) / 3.0;
			knotVector[degree + 2 * i + 1] = (uk[i] + 2 * uk[i + 1]) / 3.0;
		}
		knotVector[4] = uk[1] / 2.0;
		knotVector[knotVector.size() - degree - 2] = (uk[size - 1] + 1.0) / 2.0;
		break;
	}
	default:
	{
		std::vector<double> uk2(2 * size);
		for (int i = 0; i < size - 1; i++)
		{
			uk2[2 * i] = uk[i];
			uk2[2 * i + 1] = (uk[i] + uk[i + 1]) / 2.0;
		}
		uk2[uk2.size() - 2] = (uk2[uk2.size() - 1] + uk2[uk2.size() - 3]) / 2.0;
		knotVector = Interpolation::AverageKnotVector(degree, uk2);
	}
	}

	std::vector<std::vector<double>> A(n, std::vector<double>(n));
	for (int i = 1; i < size - 1; i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, uk[i]);
		std::vector<double> basis = Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i]);
		std::vector<std::vector<double>> derBasis = Polynomials::BasisFunctionsDerivatives(spanIndex, degree, 1, knotVector, uk[i]);
		for (int j = 0; j <= degree; j++)
		{
			A[2 * i][spanIndex - degree + j] = basis[j];
			A[2 * i + 1][spanIndex - degree + j] = derBasis[1][j];
		}
	}
	A[0][0] = 1.0;
	A[1][0] = -1.0;
	A[1][1] = 1.0;
	A[A.size() - 2][A[0].size() - 2] = -1.0;
	A[A.size() - 2][A[0].size() - 1] = 1.0;
	A[A.size() - 1][A[0].size() - 1] = 1.0;

	std::vector<std::vector<double>> right(n, std::vector<double>(3));
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			right[2 * i][j] = throughPoints[i][j];
			right[2 * i + 1][j] = unitTangents[i][j] * d;
		}
	}

	double d0 = knotVector[degree + 1] / degree;
	double dn = (1 - knotVector[knotVector.size() - degree - 2]) / degree;

	XYZ dp0 = unitTangents[0];
	XYZ dpn = unitTangents[unitTangents.size() - 1];
	XYZ qpn = throughPoints[size - 1];
	for (int j = 0; j < 3; j++)
	{
		right[1][j] = d0 * dp0[j] * d;
		right[A[0].size() - 2][j] = dn * dpn[j] * d;
		right[A[0].size() - 1][j] = qpn[j];
	}

	std::vector<std::vector<double>> result = MathUtils::SolveLinearSystem(A, right);
	for (int i = 0; i < result.size(); i++)
	{
		XYZ temp = XYZ(0, 0, 0);
		for (int j = 0; j < 3; j++)
		{
			temp[j] = result[i][j];
		}
		controlPoints[i] = XYZW(temp, 1.0);
	}
}

bool LNLib::NurbsCurve::CubicLocalInterpolation(const std::vector<XYZ>& throughPoints, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	int degree = 3;

	std::vector<XYZ> tangents;
	bool hasTangents = Interpolation::ComputeTangent(throughPoints, tangents);
	if (!hasTangents) return false;

	int size = throughPoints.size();
	std::vector<double> uk(size);
	int n = size - 1;
	uk[0] = 0;

	std::vector<XYZW> tempControlPoints;
	for (int k = 0; k < n; k++)
	{
		XYZ t0 = tangents[k];
		XYZ t3 = tangents[k + 1];
		XYZ p0 = throughPoints[k];
		XYZ p3 = throughPoints[k + 1];

		double a = 16 - (t0 + t3).SqrLength();
		double b = 12 * (p3 - p0).DotProduct(t0 + t3);
		double c = -36 * (p3 - p0).SqrLength();

		double alpha = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);

		XYZ pk0 = throughPoints[k];
		XYZ pk1 = p0 + (1 / 3) * alpha * t0;
		XYZ pk2 = p3 - (1 / 3) * alpha * t3;

		uk[k + 1] = uk[k] + 3 * (pk1 - pk0).Length();

		tempControlPoints.emplace_back(XYZW(pk1, 1));
		tempControlPoints.emplace_back(XYZW(pk2, 1));
	}

	int kvSize = 2 * (degree + 1) + 2 * (size - 1);
	knotVector.resize(kvSize);
	for (int i = 0; i <= degree; i++)
	{
		knotVector[i] = 0;
		knotVector[kvSize - 1 - i] = 1;
	}

	for (int i = 1; i <= n - 1; i = i + 2)
	{
		knotVector[degree + i] = knotVector[degree + (i + 1)] = uk[i] / uk[n];
	}

	int tSize = tempControlPoints.size();
	controlPoints.resize(tSize + 2);
	controlPoints[0] = XYZW(throughPoints[0], 1);
	controlPoints[tSize + 2 - 1] = XYZW(throughPoints[n], 1);
	for (int i = 0; i < tSize; i++)
	{
		controlPoints[i + 1] = tempControlPoints[i];
	}

	return true;
}

void LNLib::NurbsCurve::LeastSquaresApproximation(int degree, const std::vector<XYZ>& throughPoints, int controlPointsCount, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(controlPointsCount > 0, "controlPointsCount", "controlPointsCount must greater than zero.");
	int size = throughPoints.size();
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, size, controlPointsCount), "controlPointsCount", "Arguments must fit: m = n + p + 1");

	int m = size - 1;
	int n = controlPointsCount - 1;
	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	knotVector = Interpolation::ComputeKnotVector(degree, size, controlPointsCount, uk);

	std::vector<XYZ> Rk(m);
	Rk[0] = XYZ(0, 0, 0);
	for (int k = 1; k <= m - 1; k++)
	{
		double N0p = Polynomials::OneBasisFunction(0, degree, knotVector, uk[k]);
		double Nnp = Polynomials::OneBasisFunction(n, degree, knotVector, uk[k]);
		Rk[k] = throughPoints[k] - N0p * throughPoints[0] - Nnp * throughPoints[m];
	}

	std::vector<std::vector<double>> N(m - 1, std::vector<double>(n - 1));
	for (int i = 0; i <= m - 2; i++)
	{
		for (int j = 0; j <= n - 2; j++)
		{
			N[i][j] = Polynomials::OneBasisFunction(j + 1, degree, knotVector, uk[i + 1]);
		}
	}
	std::vector<std::vector<double>> Nt;
	MathUtils::Transpose(N, Nt);
	std::vector<std::vector<double>> A = MathUtils::MatrixMultiply(Nt, N);

	std::vector<XYZ> R(n - 1);
	for (int i = 0; i <= n - 2; i++)
	{
		XYZ temp = XYZ(0, 0, 0);
		for (int j = 1; j <= m - 1; j++)
		{
			temp += Polynomials::OneBasisFunction(i + 1, degree, knotVector, uk[j]) * Rk[j];
		}
		R[i] = temp;
	}

	/*std::vector<XYZ> tempControlPoints = Interpolation::ComputerControlPointsByLUDecomposition(A, R);
	for (int i = 0; i < size; i++)
	{
		controlPoints[i] = XYZW(tempControlPoints[i], 1);
	}*/
}

bool LNLib::NurbsCurve::WeightedAndContrainedLeastSquaresApproximation(unsigned int degree, const std::vector<XYZ>& throughPoints, const std::vector<double>& weights, const std::vector<XYZ>& tangents, const std::vector<int>& tangentIndices, const std::vector<double>& weightedTangents, int controlPointsCount, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	int n = controlPointsCount - 1;
	int ru = -1;
	int rc = -1;
	int size = static_cast<int>(throughPoints.size());
	int r = size - 1;
	for (int i = 0; i <= r; i++)
	{
		if (MathUtils::IsGreaterThan(weights[i], 0.0))
		{
			ru = ru + 1;
		}
		else
		{
			rc = rc + 1;
		}
	}
	int su = -1;
	int sc = -1;

	int dsize = static_cast<int>(tangents.size());
	int s = dsize - 1;
	for (int j = 0; j <= s; j++)
	{
		if (MathUtils::IsGreaterThan(weights[j], 0.0))
		{
			su = su + 1;
		}
		else
		{
			sc = sc + 1;
		}
	}
	int mu = ru + su + 1;
	int mc = rc + sc + 1;
	if (mc >= n || mc + n >= mu + 1)
	{
		return false;
	}

	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	knotVector = Interpolation::ComputeKnotVector(degree, size, controlPointsCount, uk);

	int j = 0;
	int mu2 = 0;
	int mc2 = 0;

	std::vector<std::vector<double>> N(mu + 1);
	std::vector<XYZ> S(mu + 1);
	std::vector<double> W(mu + 1);
	std::vector<std::vector<double>> M(mc + 1);
	std::vector<XYZ> T(mc + 1);
	std::vector<XYZ> A(mc + 1);

	for (int i = 0; i <= r; i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, uk[i]);
		int dflag = 0;
		if (j <= s)
		{
			if (i == tangentIndices[j])
			{
				dflag = 1;
			}
		}
		std::vector<std::vector<double>> basis;
		if (dflag == 0)
		{
			basis[i] = Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i]);
		}
		else
		{
			basis = Polynomials::BasisFunctionsDerivatives(spanIndex, degree, 1, knotVector, uk[i]);
		}
		if (MathUtils::IsGreaterThan(weights[i], 0.0))
		{
			W[mu2] = weights[i];
			N[mu2] = basis[0];
			S[mu2] = W[mu2] * throughPoints[i];
			mu2 = mu2 + 1;
		}
		else
		{
			basis[0] = M[mc2];
			T[mc2] = throughPoints[i];
			mc2 = mc2 + 1;
		}
		if (dflag == 1)
		{
			if (MathUtils::IsGreaterThan(weightedTangents[j], 0.0))
			{
				N[mu2] = basis[1];
				S[mu2] = W[mu2] * tangents[j];
				mu2 = mu2 + 1;
			}
			else
			{
				M[mc2] = basis[1];
				T[mc2] = tangents[j];
				mc2 = mc2 + 1;
			}
			j = j + 1;
		}
	}
	std::vector<std::vector<double>> NT;
	MathUtils::Transpose(N, NT);

	std::vector<std::vector<double>> NTW = MathUtils::MatrixMultiply(NT, MathUtils::MakeDiagonal(W.size()));
	std::vector<std::vector<double>> NTWN = MathUtils::MatrixMultiply(NTW, N);

	std::vector<XYZ> NTWS = Interpolation::ComputerMatrixMultiplyPoints(NTW, S);

	/*std::vector<XYZ> tempControlPoints;
	if (mc < 0)
	{
		tempControlPoints = Interpolation::ComputerControlPointsByLUDecomposition(NTWN, NTWS);
	}
	else
	{
		std::vector<std::vector<double>> inverseNTWN;
		if (!MathUtils::MakeInverse(NTWN, inverseNTWN))
		{
			return false;
		}
		std::vector<std::vector<double>> MT;
		MathUtils::Transpose(M, MT);
		std::vector<std::vector<double>> MinverseNTWN = MathUtils::MatrixMultiply(M, inverseNTWN);
		std::vector<std::vector<double>> MinverseNTWNMT = MathUtils::MatrixMultiply(MinverseNTWN, MT);

		std::vector<XYZ> MinverseNTWNNTWS = Interpolation::ComputerMatrixMultiplyPoints(MinverseNTWN, NTWS);
		std::vector<XYZ> MinverseNTWNNTWSminusT(mc + 1);
		for (int i = 0; i <= mc; i++)
		{
			MinverseNTWNNTWSminusT.emplace_back(MinverseNTWNNTWS[i] - T[i]);
		}

		std::vector<XYZ> A = Interpolation::ComputerControlPointsByLUDecomposition(MinverseNTWNMT, MinverseNTWNNTWSminusT);
		std::vector<XYZ> MTA = Interpolation::ComputerMatrixMultiplyPoints(MT, A);
		std::vector<XYZ> NTWSminusMTA(n + 1);
		for (int i = 0; i <= n; i++)
		{
			NTWSminusMTA.emplace_back(NTWS[i] - MTA[i]);
		}
		tempControlPoints = Interpolation::ComputerMatrixMultiplyPoints(inverseNTWN, NTWSminusMTA);
	}

	for (int i = 0; i < static_cast<int>(tempControlPoints.size()); i++)
	{
		controlPoints[i] = XYZW(tempControlPoints[i], 1);
	}*/
	return true;
}

double LNLib::NurbsCurve::ComputerRemoveKnotErrorBound(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, int removalIndex)
{
	int ord = static_cast<int>(degree + 1);
	int r = removalIndex;
	double u = knotVector[r];
	int s = Polynomials::GetKnotMultiplicity(knotVector, u);
	int last = r - s;
	int first = static_cast<int>(r - degree);
	int off = first - 1;

	std::vector<XYZW> temp;
	temp.resize(last + 1 - off + 1);
	temp[0] = controlPoints[off];
	temp[last + 1 - off] = controlPoints[last + 1];

	int i = first, j = last;
	int ii = 1, jj = last - off;

	double alfi = 0.0, alfj = 0.0;

	while (j - i > 0)
	{
		alfi = (u - knotVector[i]) / (knotVector[i + ord] - knotVector[i]);
		alfj = (u - knotVector[j]) / (knotVector[j + ord] - knotVector[j]);
		temp[ii] = (controlPoints[i] - (1.0 - alfi) * temp[ii - 1]) / alfi;
		temp[jj] = (controlPoints[j] - alfj * temp[jj + 1]) / (1.0 - alfj);

		i += 1;
		ii += 1;
		j = j - 1;
		jj = jj - 1;
	}
	if (j - i < 0)
	{
		return temp[ii - 1].Distance(temp[jj + 1]);
	}
	else
	{
		alfi = (u - knotVector[i]) / (knotVector[i + ord] - knotVector[i]);
		return controlPoints[i].Distance(alfi * temp[ii + 1] + (1.0 - alfi) * temp[ii - 1]);
	}
}

void LNLib::NurbsCurve::RemoveKnotsByGivenBound(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, const std::vector<double> params, std::vector<double>& error, double maxError, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	int knotSize = static_cast<int>(knotVector.size());
	std::vector<double> Br(knotSize, Constants::MaxDistance);
	std::vector<int> S(knotSize, 0);
	std::vector<int> Nl(knotSize);
	std::vector<int> Nr(knotSize, params.size() - 1);

	std::vector<double> uk = params;
	int ukSize = static_cast<int>(uk.size());
	std::vector<double> NewError(ukSize);
	std::vector<double> temp(ukSize);

	int i, BrMinIndex = 0;
	int r, s, Rstart, Rend, k;
	double BrMin, u;

	int controlPointsSize = static_cast<int>(controlPoints.size());
	int n = controlPointsSize - 1;
	s = 1;
	for (int i = degree + 1; i < controlPointsSize; i++)
	{
		if (knotVector[i] < knotVector[i + 1])
		{
			Br[i] = ComputerRemoveKnotErrorBound(degree, knotVector, controlPoints, i);
			S[i] = Polynomials::GetKnotMultiplicity(knotVector, knotVector[i]);
			s = 1;
		}
		else {
			Br[i] = Constants::MaxDistance;
			S[i] = 1;
			s++;
		}
	}

	Nl[0] = 0;
	for (int i = 0; i < ukSize; i++)
	{
		int np = static_cast<int>(knotVector.size() - degree) - 2;
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, uk[i]);
		if (!Nl[spanIndex])
		{
			Nl[spanIndex] = i;
		}
		if (i + 1 < ukSize)
		{
			Nr[spanIndex] = i + 1;
		}
	}

	std::vector<double> tempU = knotVector;
	std::vector<XYZW> tempCP = controlPoints;

	while (1)
	{
		double minStandard = Constants::MaxDistance;
		for (int i = 0; i < Br.size(); i++)
		{
			if (Br[i] < minStandard)
			{
				BrMinIndex = i;
				minStandard = Br[i];
			}
		}
		BrMin = Br[BrMinIndex];

		if (BrMin == Constants::MaxDistance)break;
		r = BrMinIndex;
		s = S[BrMinIndex];

		Rstart = std::max(r - degree, degree + 1);
		Rend = std::min(r + static_cast<int>(degree) - S[r + degree] + 1, n);
		Rstart = Nl[Rstart];
		Rend = Nr[Rend];

		bool removable = true;
		for (i = Rstart; i <= Rend; i++)
		{
			double a;
			s = S[r];
			if ((degree + s) % 2)
			{
				u = uk[i];
				k = (degree + s + 1) / 2;
				a = tempU[r] - tempU[r - k + 1];
				a /= tempU[r - k + degree + 2] - tempU[r - k + 1];
				NewError[i] = (1.0 - a) * Br[r] * Polynomials::OneBasisFunction(r - k + 1, degree, tempU, u);
			}
			else
			{
				u = uk[i];
				k = (degree + s) / 2;
				NewError[i] = Br[r] * Polynomials::OneBasisFunction(r - k, degree, tempU, u);
			}
			temp[i] = NewError[i] + error[i];
			if (temp[i] > maxError)
			{
				removable = false;
				Br[r] = Constants::MaxDistance;
				break;
			}
		}

		if (removable)
		{
			std::vector<double> tempNewU;
			std::vector<XYZW> tempNewCP;
			RemoveKnot(degree, tempU, tempCP, r, 1, tempNewU, tempNewCP);
			controlPointsSize = tempNewCP.size();
			for (int i = Rstart; i <= Rend; i++)
			{
				error[i] = temp[i];
			}

			if (controlPointsSize <= degree + 1)
			{
				break;
			}

			Rstart = Nl[r - degree - 1];
			Rend = Nr[r - S[r]];
			int spanIndex = 0;
			int oldspanIndex = -1;
			for (int k = Rstart; k <= Rend; k++)
			{
				int np = static_cast<int>(tempNewU.size() - degree) - 2;
				spanIndex = Polynomials::GetKnotSpanIndex(degree, tempNewU, uk[i]);
				if (spanIndex != oldspanIndex)
				{
					Nl[spanIndex] = k;
				}
				if (k + 1 < ukSize)
				{
					Nr[spanIndex] = k + 1;
				}
				oldspanIndex = spanIndex;
			}
			for (int k = r - S[r] + 1; k < Nl.size() - 1; k++) {
				Nl[k] = Nl[k + 1];
				Nr[k] = Nr[k + 1];
			}
			Nl.resize(Nl.size() - 1);
			Nr.resize(Nr.size() - 1);

			Rstart = std::max(r - degree, degree + 1);
			Rend = std::min(r + static_cast<int>(degree) - S[r] + 1, controlPointsSize);
			s = S[Rstart];
			for (i = Rstart; i <= Rend; i++)
			{
				if (tempNewU[i] < tempNewU[i + 1])
				{
					Br[i] = ComputerRemoveKnotErrorBound(degree, tempNewU, tempNewCP, i);
					S[i] = s;
					s = 1;
				}
				else {
					Br[i] = Constants::MaxDistance;
					S[i] = 1;
					s++;
				}
			}
			for (int i = Rend + 1; i < static_cast<int>(Br.size() - 1); i++)
			{
				Br[i] = Br[i + 1];
				S[i] = S[i + 1];
			}
			Br.resize(Br.size() - 1);

			tempU = tempNewU;
			tempCP = tempNewCP;
		}
		else
		{
			Br[r] = Constants::MaxDistance;
		}
	}
}

void LNLib::NurbsCurve::GlobalCurveApproximationByErrorBound(unsigned int degree, const std::vector<XYZ>& throughPoints, double maxError, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	int size = static_cast<int>(throughPoints.size());
	int m = size - 1;
	std::vector<double> error(size, 0.0);

	std::vector<double> U;
	U.resize(2 + m - 1 + 2);
	std::vector<XYZW> P;
	P.resize(size);

	int tempDeg = 1;
	for (int i = 0; i < uk.size(); i++)
	{
		U[i + tempDeg] = uk[i];
	}
	U[0] = 0;
	U[U.size() - 1] = 1.0;

	for (int i = 0; i < size; i++)
	{
		P[i] = XYZW(throughPoints[i], 1);
	}

	int n = m;
	int nh = -1;
	std::vector<double> Uh;
	std::vector<XYZW> Ph;

	for (int deg = 1; deg <= degree + 1; deg++)
	{
		RemoveKnotsByGivenBound(deg, U, P, uk, error, maxError, Uh, Ph);
		nh = static_cast<int>(Ph.size() - 1);

		if (deg == degree) break;

		std::vector<double> tU;
		std::vector<XYZW> tP;
		ElevateDegree(deg, Uh, Ph, 1, tU, tP);
		U = tU;
		P = tP;
		n = static_cast<int>(P.size() - 1);

		LeastSquaresApproximation(deg + 1, throughPoints, P.size(), tU, tP);
		U = tU;
		P = tP;

		for (int i = 0; i < size; i++)
		{
			double param = GetParamOnCurve(deg + 1, U, P, throughPoints[i]);
			uk[i] = param;
			XYZ p = GetPointOnCurve(deg + 1, U, param, P);
			error[i] = p.Distance(throughPoints[i]);
		}
	}

	if (n == nh)
	{
		knotVector = U;
		controlPoints = P;
		return;
	}
	else
	{
		std::vector<double> tU;
		std::vector<XYZW> tP;
		LeastSquaresApproximation(degree, throughPoints, P.size(), tU, tP);
		for (int i = 0; i < size; i++)
		{
			double param = GetParamOnCurve(degree, U, P, throughPoints[i]);
			uk[i] = param;
			XYZ p = GetPointOnCurve(degree, U, param, P);
			error[i] = p.Distance(throughPoints[i]);
		}

		RemoveKnotsByGivenBound(degree, tU, tP, uk, error, maxError, Uh, Ph);
		knotVector = Uh;
		controlPoints = Ph;
	}
}

bool LNLib::NurbsCurve::LocalRationalQuadraticCurveApproximation(int startPointIndex, int endPointIndex, const std::vector<XYZ>& throughPoints, const XYZ& startTangent, const XYZ& endTangent, double maxError, std::vector<XYZW>& middleControlPoints)
{
	XYZ startPoint = throughPoints[startPointIndex];
	XYZ endPoint = throughPoints[endPointIndex];
	if (endPointIndex - startPointIndex == 1)
	{
		return BezierCurve::ComputerMiddleControlPointsOnQuadraticCurve(startPoint, startTangent, endPoint, endTangent, middleControlPoints);
	}
	double alf1, alf2 = 0.0;
	XYZ R(0, 0, 0);
	CurveCurveIntersectionType type = Intersection::ComputeRays(startPoint, startTangent, endPoint, endTangent, alf1, alf2, R);
	if (type == CurveCurveIntersectionType::Coincident)
	{
		middleControlPoints.emplace_back(XYZW((startPoint + endPoint) / 2, 1));
		return true;
	}
	else if (type == CurveCurveIntersectionType::Skew ||
		type == CurveCurveIntersectionType::Parallel)
	{
		return false;
	}
	if (MathUtils::IsLessThanOrEqual(alf1, 0.0) ||
		MathUtils::IsGreaterThanOrEqual(alf2, 0.0))
	{
		return false;
	}
	double s = 0.0;
	XYZ V = endPoint - startPoint;
	for (int i = startPointIndex + 1; i <= endPointIndex - 1; i++)
	{
		XYZ V1 = throughPoints[i] - R;
		type = Intersection::ComputeRays(startPoint, V, R, V1, alf1, alf2, R);
		if (type == CurveCurveIntersectionType::Intersecting)
		{
			if (MathUtils::IsLessThanOrEqual(alf1, 0.0) ||
				MathUtils::IsGreaterThanOrEqual(alf1, 1.0) ||
				MathUtils::IsLessThanOrEqual(alf2, 0.0))
			{
				return false;
			}
			XYZ S = (1 - s) * (startPoint + endPoint) / 2 + s * R;
			double wi = 0.0;
			if (CreateOneConicArc(startPoint, V, R, V1, S, R, wi))
			{
				s = s + wi / (1 + wi);
			}
			else
			{
				return false;
			}
		}
		s = s / (endPointIndex - startPointIndex - 1);
		double w = s / (1.0 - s);

		std::vector<XYZW> controlPoints = { XYZW(startPoint,1), XYZW(R,w), XYZW(endPoint,1) };
		for (int i = startPointIndex + 1; i < endPointIndex - 1; i++)
		{
			XYZ tp = throughPoints[i];
			double minValue = Constants::MaxDistance;
			for (int i = 0; i < 99; i++)
			{
				double currentU = 0.01 * i;
				XYZ currentPoint;
				currentPoint = BezierCurve::GetPointOnQuadraticArc(controlPoints[0], controlPoints[1], controlPoints[2], currentU);

				double nextU = 0.01 * (i + 1);
				XYZ nextPoint;
				nextPoint = BezierCurve::GetPointOnQuadraticArc(controlPoints[0], controlPoints[1], controlPoints[2], nextU);

				XYZ vector1 = currentPoint - tp;
				XYZ vector2 = nextPoint - currentPoint;
				double dot = vector1.DotProduct(vector2);

				XYZ projectPoint;
				if (dot < 0)
				{
					projectPoint = currentPoint;
				}
				else if (dot > 1)
				{
					projectPoint = nextPoint;
				}
				else
				{
					projectPoint = currentPoint + dot * vector1.Normalize();
				}
				double distance = (tp - projectPoint).Length();
				if (distance < minValue)
				{
					minValue = distance;
				}
			}
			if (MathUtils::IsGreaterThan(minValue, maxError))
			{
				return false;
			}
		}
		middleControlPoints.emplace_back(XYZW(R, w));
		return true;
	}
	return false;
}

bool LNLib::NurbsCurve::LocalNonRationalCubicCurveApproximation(int startPointIndex, int endPointIndex, const std::vector<XYZ>& throughPoints, const XYZ& startTangent, const XYZ& endTangent, double maxError, std::vector<XYZW>& middleControlPoints)
{
	int size = static_cast<int>(throughPoints.size());
	XYZ startPoint = throughPoints[startPointIndex];
	XYZ endPoint = throughPoints[endPointIndex];
	double alpha = 0.0;
	double beta = 0.0;
	if (endPointIndex - startPointIndex == 1)
	{
		XYZ Dks(0, 0, 0);
		XYZ Dke(0, 0, 0);
		if (startPointIndex == 0)
		{
			Dks = startTangent;
		}
		if (endPointIndex == size - 1)
		{
			Dke = endTangent;
		}
		if (startPointIndex != 0 &&
			endPointIndex != size - 1)
		{
			Dks = endPoint.Distance(startPoint) / startPoint.Distance(throughPoints[startPointIndex - 1]) * startTangent;
			Dks = throughPoints[endPointIndex + 1].Distance(endPoint) / endPoint.Distance(startPoint) * endTangent;
		}
		alpha = Dks.Length() / 3;
		beta = -Dke.Length() / 3;
		XYZW P1((startPoint + alpha * startTangent), 1);
		XYZW P2((endPoint + beta * endTangent), 1);
		middleControlPoints.emplace_back(P1);
		middleControlPoints.emplace_back(P2);
		return true;
	}
	int dk = endPointIndex - startPointIndex;
	bool isLine = true;
	XYZ tempStandard(0, 0, 0);
	for (int i = startPointIndex; i < dk + 1; i++)
	{
		XYZ tp = throughPoints[i];
		XYZ direction = (tp - startPoint).Normalize();
		if (tempStandard.IsZero())
		{
			tempStandard = direction;
		}
		else
		{
			if (!tempStandard.IsAlmostEqualTo(direction))
			{
				isLine = false;
				break;
			}
		}
	}
	if (isLine)
	{
		XYZW P1((2 * startPoint + endPoint) / 3, 1);
		XYZW P2((startPoint + 2 * endPoint) / 3, 1);
		middleControlPoints.emplace_back(P1);
		middleControlPoints.emplace_back(P2);
		return true;
	}
	std::vector<double> uh;
	std::vector<double> alfak;
	std::vector<double> betak;
	for (int k = 1; k < dk; k++)
	{
		XYZ normalPi = (endPoint - startPoint).Normalize().CrossProduct(startTangent);
		XYZ t = throughPoints[k + startPointIndex];
		XYZ tt = (t - startPoint).Normalize().CrossProduct(endTangent);
		if (normalPi.Normalize().IsAlmostEqualTo(tt.Normalize()) ||
			normalPi.Normalize().IsAlmostEqualTo(tt.Normalize().Negative()))
		{
			uh = Interpolation::GetChordParameterization(throughPoints, startPointIndex, endPointIndex);
			double ak = 0.0;
			double bk = 0.0;
			double s = 1 - uh[k];
			double t = uh[k];

			XYZ a1 = 3 * pow(s, 2) * t * startTangent;
			XYZ b1 = 3 * s * pow(t, 2) * endTangent;
			XYZ c1 = throughPoints[k] - (pow(s, 3) + 3 * pow(s, 2) * t) * startPoint - (pow(t, 3) + 3 * s * pow(t, 2)) * endPoint;

			double alk = (uh[k] - uh[k - 1]) / (uh[k + 1] - uh[k - 1]);
			XYZ dk = (throughPoints[k] - throughPoints[k - 1]) / (uh[k] - uh[k - 1]);
			XYZ dk1 = (throughPoints[k + 1] - throughPoints[k]) / (uh[k + 1] - uh[k]);
			XYZ Dk = ((1 - alk) * dk + alk * dk1).Normalize();

			XYZ a2 = s * (s - 2 * t) * (Dk.CrossProduct(startTangent));
			XYZ b2 = t * (2 * s - t) * (Dk.CrossProduct(endTangent));
			XYZ c2 = 2 * s * t * (Dk.CrossProduct(startPoint - endPoint));

			ak = (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
			bk = (c1 * a2 - c2 * a1) / (b1 * a2 - b2 * a1);

			if (MathUtils::IsGreaterThan(ak, 0.0) &&
				MathUtils::IsLessThan(bk, 0.0))
			{
				alfak[k] = ak;
				betak[k] = bk;
			}
			else
			{
				return false;
			}
		}
		else
		{
			XYZ Pd(0, 0, 0);
			LinePlaneIntersectionType type0 = Intersection::ComputeLineAndPlane(normalPi, startPoint, throughPoints[k + startPointIndex], startTangent, Pd);
			if (!(type0 == LinePlaneIntersectionType::Intersecting)) return false;
			double param0 = 0.0; double param1 = 0.0;
			XYZ Pc(0, 0, 0);
			CurveCurveIntersectionType type1 = Intersection::ComputeRays(startPoint, endPoint - startPoint, Pd, startTangent, param0, param1, Pc);
			if (!(type1 == CurveCurveIntersectionType::Intersecting)) return false;
			double gamma = Pc.Distance(endPoint) / startPoint.Distance(endPoint);
			if (MathUtils::IsLessThan(gamma, 0.0) ||
				MathUtils::IsGreaterThan(gamma, 1.0))
			{
				return false;
			}
			uh[k] = MathUtils::ComputerCubicEquationsWithOneVariable(2, -3, 0, 1 - gamma);
			if (MathUtils::IsLessThan(uh[k], 0.0) ||
				MathUtils::IsGreaterThan(uh[k], 1.0))
			{
				return false;
			}
			else
			{
				double a = Pc.Distance(Pd);
				double b = -Pd.Distance(throughPoints[k + startPointIndex]);
				double b13 = Polynomials::Bernstein(1, 3, uh[k]);
				double b23 = Polynomials::Bernstein(2, 3, uh[k]);
				alfak[k] = a / b13;
				betak[k] = b / b23;
			}
		}
	}
	alpha = 0.0;
	beta = 0.0;
	for (int k = 1; k < dk; k++)
	{
		alpha += alfak[k];
		beta += betak[k];
	}
	alpha = alpha / (dk - 1);
	beta = beta / (dk - 1);
	XYZ P1 = throughPoints[startPointIndex] + alpha * startTangent;
	XYZ P2 = throughPoints[endPointIndex] + beta * endTangent;
	middleControlPoints.emplace_back(XYZW(P1, 1));
	middleControlPoints.emplace_back(XYZW(P2, 1));
	for (int k = 1; k < dk; k++)
	{
		double u = uh[k];
		double deltaAk = alfak[k] - alpha;
		double deltaBk = betak[k] - beta;
		double e = (deltaAk * Polynomials::Bernstein(1, 3, u) * startTangent).Distance(deltaBk * Polynomials::Bernstein(2, 3, u) * endTangent);
		if (MathUtils::IsLessThanOrEqual(e, maxError)) continue;
		else
		{
			std::vector<XYZW> cps = { XYZW(startPoint,1), XYZW(P1,1), XYZW(P2,1),  XYZW(endPoint,1) };
			double tempParam = GetParamOnCurve(3, uh, cps, throughPoints[startPointIndex + k]);
			std::vector<XYZ> cpts = { startPoint, P1, P2, endPoint };
			XYZ t = BezierCurve::GetPointOnCurveByBernstein(3, cpts, tempParam);
			XYZ p = BezierCurve::GetPointOnCurveByBernstein(3, cpts, u);
			double ek = p.Distance(t);
			if (MathUtils::IsGreaterThan(ek, maxError))
			{
				return false;
			}
		}
	}
	return true;
}

void LNLib::NurbsCurve::ControlPointReposition(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, XYZW newControlPoint, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	int size = controlPoints.size();
	int n = size - 1;

	double ti = GetNode(degree, knotVector, n);

	int k = 0;
	int kvSize = knotVector.size();
	double diff = Constants::MaxDistance;
	for (int i = 0; i < kvSize; i++)
	{
		int endIndex = i + degree + 1;
		if (endIndex > kvSize)
		{
			break;
		}
		else
		{
			for (int j = 1; j <= j + degree + 1; j++)
			{
				double d = abs(knotVector[j] - ti);
				if (d < diff)
				{
					diff = d;
					k = j;
				}
			}
		}
	}

	double u = knotVector[k];
	double tk = GetNode(degree, knotVector, k);
	double tk1 = GetNode(degree, knotVector, k + 1);

	if (MathUtils::IsGreaterThanOrEqual(abs(u - tk), abs(tk1 - u)))
	{
		double sum = 0.0;
		for (int j = 2; j <= degree; j++)
		{
			sum += knotVector[k + j];
		}
		double insertKnot = degree * u - sum;
		InsertKnot(degree, knotVector, controlPoints, insertKnot, 1, updatedKnotVector, updatedControlPoints);
		updatedControlPoints[k + 1] = newControlPoint;
	}
	else
	{
		updatedKnotVector = knotVector;
		updatedControlPoints = controlPoints;
		updatedControlPoints[k] = newControlPoint;
	}
}

double LNLib::NurbsCurve::WeightModification(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, XYZ pointOnCurve, XYZ selectedControlPoint, XYZ newControlPoint, bool isPull)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	XYZ dir1 = (selectedControlPoint - newControlPoint).Normalize();
	XYZ dir2 = (newControlPoint - pointOnCurve).Normalize();
	if (!dir1.IsZero() && !dir2.IsZero())
	{
		VALIDATE_ARGUMENT(dir1.IsAlmostEqualTo(dir2) || dir1.IsAlmostEqualTo(dir2.Negative()), "newControlPoint", "NewControlPoint must on the line created by pointOnCurve and selectedControlPoint");
	}

	double u = GetParamOnCurve(degree, knotVector, controlPoints, pointOnCurve);
	double distance = pointOnCurve.Distance(newControlPoint);
	double pkp = selectedControlPoint.Distance(pointOnCurve);
	int k = Polynomials::GetKnotSpanIndex(degree, knotVector, u);
	double wk = controlPoints[k].GetW();

	double nip = Polynomials::OneBasisFunction(k, degree, knotVector, u) * wk;
	double sum = 0.0;
	for (int r = 0; r < controlPoints.size(); r++)
	{
		sum += Polynomials::OneBasisFunction(r, degree, knotVector, u) * controlPoints[r].GetW();
	}
	double rkp = nip / sum;
	return isPull ? wk * (1 + (abs(distance) / (rkp) * (pkp - distance))) : wk * (1 - (abs(distance) / (rkp) * (pkp + distance)));
}

std::vector<LNLib::XYZW> LNLib::NurbsCurve::Warping(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double warpDistance, const XYZ& planeNormal)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	// to be continue...
	return std::vector<XYZW>();
}

void LNLib::NurbsCurve::Flattening(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, XYZ lineStartPoint, XYZ lineEndPoint, double flattenStartParam, double flattenEndParam, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	VALIDATE_ARGUMENT(!lineStartPoint.IsAlmostEqualTo(lineEndPoint), "lineEndPoint", "lineEndPoint must not be equals to lineStartPoint.");
	VALIDATE_ARGUMENT_RANGE(flattenStartParam, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(flattenEndParam, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(flattenEndParam, flattenStartParam), "flattenEndParam", "flattenEndParam must be greater than flattenStartParam.");

	int spanMinIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, flattenStartParam);
	int spanMaxIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, flattenEndParam);
	
	std::unordered_map<int, XYZ> selectedControlPoints;
	for (int i = spanMinIndex - degree; i <= spanMaxIndex - degree - 1; i++)
	{
		XYZ p = const_cast<XYZW&>(controlPoints[i]).ToXYZ(true);
		selectedControlPoints.insert({ i, p });
	}
	
	int projectCount = 0;
	updatedControlPoints = controlPoints;
	XYZ lineVector = (lineEndPoint - lineStartPoint).Normalize();
	double lineDistance = lineStartPoint.Distance(lineEndPoint);
	for (auto it = selectedControlPoints.begin(); it != selectedControlPoints.end(); ++it)
	{
		int index = it->first;
		XYZ current = it->second;
		XYZ project = Projection::PointToRay(lineStartPoint, lineVector, current);
		double dis = current.Distance(project);
		if (MathUtils::IsLessThanOrEqual(dis, lineDistance))
		{
			projectCount++;
			updatedControlPoints[index] = XYZW(project,updatedControlPoints[index].GetW());
		}
	}
	if (projectCount >= degree + 1)
	{
		updatedKnotVector = knotVector;
	}
	else
	{
		double ui = knotVector[0];
		int kvSize = knotVector.size();
		for (int i = 0; i < kvSize; i++)
		{
			double current = knotVector[i];
			if (MathUtils::IsLessThan(current, flattenStartParam))
			{
				ui = current;
				continue;
			}
			else
			{
				break;
			}
		}

		double uj = knotVector[kvSize - 1];
		for (int j = kvSize - 1; j >= 0; j--)
		{
			double current = knotVector[j];
			if (MathUtils::IsGreaterThan(current, flattenEndParam))
			{
				uj = current;
				continue;
			}
			else
			{
				break;
			}
		}

		std::vector<double> insertElement = Polynomials::GetInsertedKnotElement(degree, knotVector, ui, uj);
		std::vector<XYZW> refinedControlPoints;
		RefineKnotVector(degree, knotVector, controlPoints, insertElement, updatedKnotVector, refinedControlPoints);

		selectedControlPoints.clear();
		spanMinIndex = Polynomials::GetKnotSpanIndex(degree, updatedKnotVector, flattenStartParam);
		spanMaxIndex = Polynomials::GetKnotSpanIndex(degree, updatedKnotVector, flattenEndParam);
		for (int i = spanMinIndex - degree; i <= spanMaxIndex - degree - 1; i++)
		{
			XYZ p = const_cast<XYZW&>(refinedControlPoints[i]).ToXYZ(true);
			selectedControlPoints.insert({ i, p });
		}

		for (auto it = selectedControlPoints.begin(); it != selectedControlPoints.end(); ++it)
		{
			int index = it->first;
			XYZ current = it->second;
			XYZ project = Projection::PointToRay(lineStartPoint, lineVector, current);
			double dis = current.Distance(project);
			if (MathUtils::IsLessThanOrEqual(dis, lineDistance))
			{
				refinedControlPoints[index] = XYZW(project, refinedControlPoints[index].GetW());
			}
		}
		updatedControlPoints = refinedControlPoints;
	}
}

std::vector<LNLib::XYZW> LNLib::NurbsCurve::Bending(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double bendStartParam, double bendEndParam, int bendCurveDegree, const std::vector<double>& bendCurveKnotVector, const std::vector<XYZW>& bendCurveControlPoints, const XYZ& bendCenter, double crossRatio)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	VALIDATE_ARGUMENT_RANGE(bendStartParam, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(bendEndParam, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(bendEndParam, bendStartParam), "bendEndParam", "bendEndParam must be greater than bendStartParam.");
	
	VALIDATE_ARGUMENT(bendCurveDegree > 0, "bendCurveDegree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(bendCurveKnotVector.size() > 0, "bendCurveKnotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(bendCurveKnotVector), "bendCurveKnotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(bendCurveControlPoints.size() > 0, "bendCurveControlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(bendCurveDegree, bendCurveKnotVector.size(), bendCurveControlPoints.size()), "bendCurveControlPoints", "Arguments must fit: m = n + p + 1");

	int spanMinIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, bendStartParam);
	int spanMaxIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, bendEndParam);

	std::vector<XYZW> updatedControlPoints = controlPoints;
	std::vector<double> insertElement = Polynomials::GetInsertedKnotElement(degree, knotVector, knotVector[spanMinIndex], knotVector[spanMaxIndex]);
	if (insertElement.size() > 0)
	{
		std::vector<double> updatedKnotVector;
		std::vector<XYZW> refinedControlPoints;
		RefineKnotVector(degree, knotVector, controlPoints, insertElement, updatedKnotVector, refinedControlPoints);

		spanMinIndex = Polynomials::GetKnotSpanIndex(degree, updatedKnotVector, bendStartParam);
		spanMaxIndex = Polynomials::GetKnotSpanIndex(degree, updatedKnotVector, bendEndParam);

		updatedControlPoints = refinedControlPoints;
	}
	
	std::unordered_map<int, XYZ> selectedControlPoints;
	for (int i = spanMinIndex - degree; i <= spanMaxIndex - degree - 1; i++)
	{
		XYZ p = const_cast<XYZW&>(updatedControlPoints[i]).ToXYZ(true);
		selectedControlPoints.insert({ i, p });
	}

	for (auto it = selectedControlPoints.begin(); it != selectedControlPoints.end(); ++it)
	{
		int index = it->first;
		XYZ current = it->second;
		double bendParam = GetParamOnCurve(bendCurveDegree, bendCurveKnotVector, bendCurveControlPoints, current);
		if (MathUtils::IsAlmostEqualTo(bendParam, Constants::DoubleEpsilon))
		{
			continue;
		}
		XYZ pointOnBendCurve = GetPointOnCurve(bendCurveDegree, bendCurveKnotVector, bendParam, bendCurveControlPoints);
		double si = bendCenter.Distance(pointOnBendCurve) / bendCenter.Distance(current);
		double ti = (crossRatio * si) / (1 + (crossRatio - 1) * si);
		XYZ project = bendCenter + (current - bendCenter) * ti;
		updatedControlPoints[index] = XYZW(project, updatedControlPoints[index].GetW());
	}
	return updatedControlPoints;
}

std::vector<LNLib::XYZW> LNLib::NurbsCurve::ConstraintBasedModification(int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, const std::vector<double>& constraintParams, const std::vector<XYZ>& derivativeConstraints, const std::vector<int>& appliedIndices, const std::vector<int>& appliedDegree, const std::vector<int>& fixedControlPointIndices)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	VALIDATE_ARGUMENT(constraintParams.size() > 0, "constraintParams", "ConstraintParams size must greater than zero.");
	VALIDATE_ARGUMENT(derivativeConstraints.size() > 0, "derivativeConstraints", "DerivativeConstraints size must greater than zero.");
	VALIDATE_ARGUMENT(appliedIndices.size() == derivativeConstraints.size(), "appliedIndices", "AppliedIndices size must equals to derivativeConstraints size.");
	VALIDATE_ARGUMENT(appliedDegree.size() == appliedIndices.size(), "appliedDegree", "AppliedDegree size must equals to appliedIndices size.");

	std::vector<double> ur = constraintParams;
	std::vector<XYZ> D = derivativeConstraints;
	std::vector<int> Dr = appliedIndices;
	std::vector<int> Dk = appliedDegree;

	int size = controlPoints.size();
	std::vector<std::vector<double>> B(D.size(), std::vector<double>(size));
	for (int i = 0; i < D.size(); i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, ur[Dr[i]]);
		std::vector<std::vector<double>> ders = Polynomials::BasisFunctionsDerivatives(spanIndex, degree, Dk[i], knotVector, ur[Dr[i]]);
		for (int j = 0; j <= degree; j++)
		{
			B[i][spanIndex - degree + j] = ders[Dk[i]][j];
		}
	}

	std::vector<int> remove(size, 1.0);
	std::vector<int> map(size);

	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < D.size(); i++)
		{
			if (MathUtils::IsGreaterThan(B[i][j] * B[i][j], 0.0))
			{
				remove[i] = 0;
				break;
			}
		}
	}

	for (int i = 0; i < fixedControlPointIndices.size(); i++)
	{
		remove[fixedControlPointIndices[i]] = 1;
	}

	int n = 0;
	for (int i = 0; i < size; i++)
	{
		if (!remove[i])
		{
			map[n] = i;
			n++;
		}
	}

	map.resize(n);

	std::vector<std::vector<double>> Bopt(D.size(), std::vector<double>(n));
	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < D.size(); i++)
		{
			Bopt[i][j] = B[i][map[j]];
		}
	}

	std::vector<std::vector<double>> Bt;
	MathUtils::Transpose(Bopt, Bt);
	std::vector<std::vector<double>> BBt;
	MathUtils::MakeInverse(MathUtils::MatrixMultiply(Bopt, Bt), BBt);

	auto A = MathUtils::MatrixMultiply(Bt, BBt);
	std::vector<std::vector<double>> dD(D.size(), std::vector<double>(3));
	for (int i = 0; i < D.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			dD[i][j] = D[i][j];
		}
	}

	auto dP = MathUtils::MatrixMultiply(A, dD);
	std::vector<XYZW> updatedControlPoints = controlPoints;
	for (int i = 0; i < map.size(); i++)
	{
		double weight = updatedControlPoints[map[i]].GetW();

		double x = dP[i][0];
		double y = dP[i][1];
		double z = dP[i][2];

		double wx = updatedControlPoints[map[i]].GetWX();
		double wy = updatedControlPoints[map[i]].GetWY();
		double wz = updatedControlPoints[map[i]].GetWZ();

		updatedControlPoints[map[i]] = XYZW(wx + weight * x, wy + weight * y, wz + weight * z, weight);
	}
	return updatedControlPoints;
}

void LNLib::NurbsCurve::ToClampCurve(int degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	std::vector<double> insertKnotVector;
	std::vector<XYZW> updatedControlPoints;
	InsertKnot(degree, knotVector, controlPoints, knotVector[degree], degree, insertKnotVector, updatedControlPoints);

	knotVector = insertKnotVector;
	controlPoints = updatedControlPoints;

	std::vector<double> insertKnotVectorNc;
	std::vector<XYZW> updatedControlPointsNc;
	InsertKnot(degree, knotVector, controlPoints, knotVector[controlPoints.size()], degree, insertKnotVectorNc, updatedControlPointsNc);

	knotVector.resize(insertKnotVectorNc.size() - degree - degree, 0.0);
	controlPoints.resize(knotVector.size() - degree - 1);
	for (int i = knotVector.size() - 1; i >= 0; i--)
	{
		knotVector[i] = insertKnotVectorNc[i + degree];
		if (i < controlPoints.size())
		{
			controlPoints[i] = updatedControlPointsNc[i + degree];
		}
	}
}

void LNLib::NurbsCurve::ToUnclampCurve(int degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contains one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must fit: m = n + p + 1");

	int n = controlPoints.size() - 1;
	for (int i = 0; i <= degree - 2; i++)
	{
		knotVector[degree - i - 1] = knotVector[degree - i] - (knotVector[n - i + 1] - knotVector[n - i]);
		int k = degree - 1;
		for (int j = i; j >= 0; j--)
		{
			double alfa = (knotVector[degree] - knotVector[k]) / (knotVector[degree + j + 1] - knotVector[k]);
			controlPoints[j] = (controlPoints[j] - alfa * controlPoints[j + 1]) / (1.0 - alfa);
			k = k - 1;
		}
	}
	knotVector[0] = knotVector[1] - (knotVector[n - degree + 2] - knotVector[n - degree + 1]);
	for (int i = 0; i <= degree - 2; i++)
	{
		knotVector[n + i + 2] = knotVector[n + i + 1] + (knotVector[degree + i + 1] - knotVector[degree + i]);
		for (int j = i; j >= 0; j--)
		{
			double alfa = (knotVector[n + 1] - knotVector[n - j]) / (knotVector[n - j + i + 2] - knotVector[n - j]);
			controlPoints[n - j] = (controlPoints[n - j] - (1.0 - alfa) * controlPoints[n - j - 1]) / alfa;
		}
	}
	knotVector[n + degree + 1] = knotVector[n + degree] + (knotVector[2 * degree] - knotVector[2 * degree - 1]);
}
