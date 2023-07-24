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
#include "BsplineCurve.h"
#include "Intersection.h"
#include "ValidationUtils.h"
#include "Interpolation.h"
#include <vector>
#include <algorithm>

void LNLib::NurbsCurve::GetPointOnCurve(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double paramT, XYZ& point)
{
	XYZW temp = XYZW();
	BsplineCurve::GetPointOnCurve(controlPoints, degree, paramT, knotVector, temp);
	point = temp.ToXYZ(true);
}

void LNLib::NurbsCurve::ComputeRationalCurveDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives)
{
	derivatives.resize(derivative + 1);

	std::vector<XYZW> ders;
	BsplineCurve::ComputeDerivatives(degree, knotVector, controlPoints, paramT, derivative, ders);

	std::vector<XYZ> Aders;
	Aders.resize(derivative + 1);
	for (int i = 0; i < static_cast<int>(ders.size()); i++)
	{
		Aders[i] = ders[i].ToXYZ(false);
	}
	std::vector<double> wders;
	wders.resize(derivative + 1);
	for (int i = 0; i < static_cast<int>(ders.size()); i++)
	{
		wders[i] = ders[i].GetW();
	}

	for (int k = 0; k <= static_cast<int>(derivative); k++)
	{
		XYZ v = Aders[k];
		for (int i = 1; i <= k; i++)
		{
			v = v - MathUtils::Binomial(k, i) * wders[i] * derivatives[k - i];
		}
		derivatives[k] = v/wders[0];
	}

}

void LNLib::NurbsCurve::InsertKnot(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double insertKnot, unsigned int times, std::vector<double>& insertedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	int np = static_cast<int>(controlPoints.size() - 1);
	int knotSpanIndex = Polynomials::GetKnotSpanIndex(np, degree, insertKnot, knotVector);
	int originMultiplicity = Polynomials::GetKnotMultiplicity(insertKnot, knotVector);

	if (originMultiplicity + times > degree)
	{
		times = degree - originMultiplicity;
	}

	insertedKnotVector.resize(knotVector.size() + times);
	for (int i = 0; i <= knotSpanIndex; i++)
	{
		insertedKnotVector[i] = knotVector[i];
	}
	for (int i = 1; i <= static_cast<int>(times); i++)
	{
		insertedKnotVector[knotSpanIndex + i] = insertKnot;
	}
	for (int i = knotSpanIndex + 1; i < static_cast<int>(knotVector.size()); i++)
	{
		insertedKnotVector[i+times] = knotVector[i];
	}

	updatedControlPoints.resize(controlPoints.size() + times);
	for (int i = 0; i <= knotSpanIndex - static_cast<int>(degree); i++)
	{
		updatedControlPoints[i] = controlPoints[i];
	}
	for (int i = knotSpanIndex - originMultiplicity; i < controlPoints.size(); i++)
	{
		updatedControlPoints[i+times] = controlPoints[i];
	}

	std::vector<XYZW> temp;
	temp.resize(degree - originMultiplicity + 1);
	for (int i = 0; i <= static_cast<int>(degree) - originMultiplicity; i++)
	{
		temp[i] = controlPoints[knotSpanIndex - degree + i];
	}

	int L = 0;
	for (int j = 1; j <= static_cast<int>(times); j++)
	{
		L = knotSpanIndex - degree + j;
		for (int i = 0; i <= static_cast<int>(degree) - j - originMultiplicity; i++)
		{
			double alpha = (insertKnot - knotVector[L + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[L + i]);
			temp[i] = alpha * temp[i + 1] + (1.0 - alpha) * temp[i];
		}
		updatedControlPoints[L] = temp[0];
		updatedControlPoints[knotSpanIndex + times - j - originMultiplicity] = temp[degree - j - originMultiplicity];
	}

	for (int i = L +1 ; i < knotSpanIndex - originMultiplicity; i++)
	{
		updatedControlPoints[i] = temp[i - L];
	}
}

void LNLib::NurbsCurve::GetPointOnCurveByInsertKnot(unsigned int degree,const std::vector<double>& knotVector, std::vector<XYZW>& controlPoints, double insertKnot, XYZ& point)
{
	int n = static_cast<int>(controlPoints.size() - 1);

	if (MathUtils::IsAlmostEqualTo(insertKnot, knotVector[0]))
	{
		point = controlPoints[0].ToXYZ(true);
		return;
	}
	if (MathUtils::IsAlmostEqualTo(insertKnot, knotVector[n + degree + 1]))
	{
		point = controlPoints[n].ToXYZ(true);
		return;
	}

	int knotSpanIndex = Polynomials::GetKnotSpanIndex(n, degree, insertKnot, knotVector);
	int originMultiplicity = Polynomials::GetKnotMultiplicity(insertKnot, knotVector);

	int times = degree - originMultiplicity;
	std::vector<XYZW> temp;
	temp.resize(times + 1);
	for (int i = 0; i <= times; i++)
	{
		temp[i] = controlPoints[knotSpanIndex - degree + i];
	}
	for (int j = 1; j <= times; j++)
	{
		for (int i = 0; i < times - j; i++)
		{
			double alpha = (insertKnot - knotVector[knotSpanIndex - degree + j + i]) / (knotVector[i + knotSpanIndex + 1] - knotVector[knotSpanIndex - degree + j + i]);
			temp[i] = alpha * temp[i + 1] + (1.0 - alpha) * temp[i];
		}
	}
	point = temp[0].ToXYZ(true);
}

void LNLib::NurbsCurve::RefineKnotVector(unsigned int degree, const std::vector<double>& knotVector,const std::vector<XYZW>& controlPoints, std::vector<double>& insertKnotElements, std::vector<double>& insertedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	int n = static_cast<int>(controlPoints.size() - 1);
	int m = n + degree + 1;
	int r = static_cast<int>(insertKnotElements.size() - 1);

	int a = Polynomials::GetKnotSpanIndex(n, degree, insertKnotElements[0], knotVector);
	int b = Polynomials::GetKnotSpanIndex(n, degree, insertKnotElements[r], knotVector);

	b = b + 1;

	insertedKnotVector.resize(knotVector.size() + insertKnotElements.size());
	for (int j = 0; j <= a; j++)
	{
		insertedKnotVector[j] = knotVector[j];
	}
	for (int j = b + degree; j <= m; j++)
	{
		insertedKnotVector[j + r + 1] = knotVector[j];
	}

	updatedControlPoints.resize(n + r + 1);
	for (int j = 0; j <= a - static_cast<int>(degree); j++)
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
		for (int l = 1; l <= static_cast<int>(degree); l++)
		{
			int ind = k - degree + 1;
			double alpha = insertedKnotVector[k + l] - insertKnotElements[j];
			if (MathUtils::IsAlmostEqualTo(abs(alpha), 0.0))
			{
				updatedControlPoints[ind - 1] = updatedControlPoints[ind];
			}
			else
			{
				alpha = alpha / (insertedKnotVector[k + 1] - knotVector[i - degree + 1]);
				updatedControlPoints[ind - 1] = alpha * updatedControlPoints[ind - 1] + (1.0 - alpha) * updatedControlPoints[ind];
			}
		}

		insertedKnotVector[k] = insertKnotElements[j];
		k = k - 1;
	}
}

void LNLib::NurbsCurve::ToBezierCurves(unsigned int degree, const std::vector<double>& knotVector,const std::vector<XYZW>& controlPoints, int& bezierCurvesCount, std::vector<std::vector<XYZW>>& decomposedControlPoints)
{
	int n = static_cast<int>(controlPoints.size() - 1);
	int m = n + degree + 1;

	int a = static_cast<int>(degree);
	int b = static_cast<int>(degree) + 1;

	int bezierCurves = static_cast<int>(knotVector.size() / (degree + 1)) - 1;
	decomposedControlPoints.resize(bezierCurves);
	for (int i = 0; i < bezierCurves; i++)
	{
		decomposedControlPoints[i].resize(degree + 1);
	}

	bezierCurvesCount = 0;
	for (int i = 0; i <= static_cast<int>(degree); i++)
	{
		decomposedControlPoints[bezierCurvesCount][i] = controlPoints[i];
	}

	while (b < m)
	{
		int i = b;
		while (b < m && MathUtils::IsAlmostEqualTo(knotVector[b + 1], knotVector[b]))
		{
			b++;
		}
		int multi = b - i + 1;
		if (multi < static_cast<int>(degree))
		{
			double numerator = knotVector[b] - knotVector[a];

			std::vector<double> alhpaVector;
			alhpaVector.resize(degree - multi);
			for (int j = degree; j > multi; j--)
			{
				alhpaVector[j - multi - 1] = numerator / (knotVector[a + j] - knotVector[a]);
			}

			int r = degree - multi;
			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				int s = multi + j;
				for (int k = degree; k >= s; k--)
				{
					double alpha = alhpaVector[k - s];
					decomposedControlPoints[bezierCurvesCount][k] = alpha * decomposedControlPoints[bezierCurvesCount][k] + (1.0 - alpha) * decomposedControlPoints[bezierCurvesCount][k - 1];
				}

				if (b < m)
				{
					decomposedControlPoints[bezierCurvesCount + 1][save] = decomposedControlPoints[bezierCurvesCount][degree];
				}
			}

			bezierCurvesCount += 1;
			if (b < m)
			{
				for (int i = static_cast<int>(degree) - multi; i <= static_cast<int>(degree); i++)
				{
					decomposedControlPoints[bezierCurvesCount][i] = controlPoints[b - degree + i];
				}

				a = b;
				b += 1;
			}
		}
	}
}

void LNLib::NurbsCurve::RemoveKnot(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double removeKnot, unsigned int times, std::vector<double>& restKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);
	int n = static_cast<int>(controlPoints.size() - 1);

	int order = degree + 1;
	int originMultiplicity = Polynomials::GetKnotMultiplicity(removeKnot, knotVector);
	int removeIndex = Polynomials::GetKnotSpanIndex(n, degree, removeKnot, knotVector);

	int first = removeIndex - degree;
	int last = removeIndex - originMultiplicity;
	
	if (times < 1)
	{
		restKnotVector = knotVector;
		updatedControlPoints = controlPoints;
		return;
	}

	restKnotVector.resize(knotVector.size());
	updatedControlPoints.resize(controlPoints.size());

	std::vector<XYZW> temp;
	temp.resize(times);

	int t = 0;
	for (t = 0; t < static_cast<int>(times); t++)
	{
		int off = first - 1;
		temp[0] = controlPoints[off];
		temp[last + 1 - off] = controlPoints[last + 1];

		int i = first;
		int j = last;

		int ii = 1;
		int jj = last - off;
		bool remflag = false;
		while (j - i > static_cast<int>(t))
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

		if (j - i < static_cast<int>(t))
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

		if (!remflag)
		{
			break;
		}
		else
		{
			i = first;
			j = last;

			while (j - i > t)
			{
				updatedControlPoints[i] = temp[i - off];
				updatedControlPoints[j] = temp[j - off];
				i = i + 1;
				j = j + 1;
			}
		}
		first = first - 1;
		last = last + 1;
	}

	if (t == 0)
	{
		return;
	}
	int m = n + degree + 1;
	for (int k = removeIndex + 1; k <= m; k++)
	{
		restKnotVector[k - t] = knotVector[k];
	}
	for (int k = 0; k < t; k++)
	{
		restKnotVector.pop_back();
	}

	int j = static_cast<int>((2 * removeIndex - originMultiplicity - degree) / 2);
	int i = j;
	
	for (int k = 1; k < t; k++)
	{
		if ( k%2 == 1)
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
	for (int k = 0; k < t; k++)
	{
		updatedControlPoints.pop_back();
	}
}

void LNLib::NurbsCurve::ElevateDegree(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, unsigned int times,  std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	if (degree + times <= degree)
	{
		updatedKnotVector = knotVector;
		updatedControlPoints = controlPoints;
		return;
	}

	int n = static_cast<int>(controlPoints.size() - 1);
	int m = n + degree + 1;

	int ph = degree + times;
	int ph2 = static_cast<int>(ph / 2);

	std::vector<std::vector<double>> bezalfs;
	bezalfs.resize(degree + times + 1);
	for (int i = 0; i < static_cast<int>(degree + times + 1); i++)
	{
		bezalfs[i].resize(degree + 1);
	}

	bezalfs[0][0] = bezalfs[ph][degree] = 1.0;

	for (int i = 1; i <= ph2; i++)
	{
		double inv = 1.0 / MathUtils::Binomial(ph, i);
		int mpi = std::min(static_cast<int>(degree), i);

		for (int j = std::max(0, i - static_cast<int>(times)); j <= mpi; j++)
		{
			bezalfs[i][j] = inv * MathUtils::Binomial(degree, j) * MathUtils::Binomial(times, i - j);
		}
	}

	for (int i = ph2 + 1; i <= ph - 1; i++)
	{
		int mpi = std::min(static_cast<int>(degree), i);
		for (int j = std::max(0, i - static_cast<int>(times)); j <= mpi; j++)
		{
			bezalfs[i][j] = bezalfs[ph-i][degree - j];
		}
	}

	int mh = ph;
	int kind = ph + 1;
	int r = -1;
	int a = static_cast<int>(degree);
	int b = static_cast<int>(degree) + 1;
	int cind = 1;
	double ua = knotVector[0];

	updatedControlPoints.resize(kind);
	updatedControlPoints[0] = controlPoints[0];

	updatedKnotVector.resize(ph + 1, ua);

	std::vector<XYZW> bpts;
	bpts.resize(degree + 1);
	for (int i = 0; i <= static_cast<int>(degree); i++)
	{
		bpts[i] = controlPoints[i];
	}

	std::vector<XYZW> nextbpts;
	nextbpts.resize(degree - 1);
	
	int nh = 0;
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
		r = static_cast<int>(degree) - mul;

		int lbz = oldr > 0? (oldr + 2) / 2: 1;
		int rbz = r>0? ph - (r + 1) / 2: ph;

		if (r > 0)
		{
			double numer = ub - ua;
			std::vector<double> alfs;
			alfs.resize(degree - 1);
			for (int k = static_cast<int>(degree); k > mul; k--)
			{
				alfs[k - mul - 1] = numer / (knotVector[a + k] - ua);
			}
			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				int s = mul + j;

				for (int k = static_cast<int>(degree); k >= s; k--)
				{
					bpts[k] = alfs[k - s] * bpts[k] + (1.0 - alfs[k - s]) * bpts[k - 1];
				}

				nextbpts[save] = bpts[degree];
			}
		}

		std::vector<XYZW> ebpts;
		ebpts.resize(degree + times + 1);
		for (int i = lbz; i <= ph; i++)
		{
			ebpts[i] = XYZW(0.0,0.0,0.0,0.0);
			int mpi = std::min(static_cast<int>(degree), i);
			for (int j = std::max(0, static_cast<int>(i - times)); j <= mpi; j++)
			{
				ebpts[i] += bezalfs[i][j] * bpts[j];
			}
		}

		if (oldr > 1)
		{
			int first = kind - 2;
			int last = kind;
			double den = ub - ua;
			double bet = (ub - updatedKnotVector[kind - 1] / den);
			
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
							double gam = (ub - updatedKnotVector[j - tr] / den);
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

				first = first - 1;
				last += 1;
			}
		}

		if (a != degree)
		{
			for (int i = 0; i < ph - oldr; i++)
			{
				updatedKnotVector[kind] = ua;
				kind += 1;
			}
		}

		for (int j = lbz; j <= rbz; j++)
		{
			updatedControlPoints[cind] = ebpts[j];
			cind += 1;
		}

		if (b < m)
		{
			for (int j = 0; j < r; j++)
			{
				bpts[j] = nextbpts[j];
			}
			for (int j = r; j <= static_cast<int>(degree); j++)
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
	nh = mh - ph - 1;
}

bool LNLib::NurbsCurve::ReduceDegree(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	bool canReduce = ValidationUtils::IsValidDegreeReduction(degree);
	if (!canReduce) return false;

	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);

	int ph = static_cast<int>(degree) - 1;
	int mh = ph;

	int kind = ph + 1;
	int r = -1;
	int a = static_cast<int>(degree);

	int b = static_cast<int>(degree) + 1;
	int cind = 1;
	int mult = static_cast<int>(degree);

	int n = static_cast<int>(controlPoints.size() - 1);
	int m = n + static_cast<int>(degree) + 1;

	std::vector<XYZW> bpts;
	bpts.resize(degree + 1);

	std::vector<XYZW> nextbpts;
	nextbpts.resize(degree - 1);

	std::vector<XYZW> rbpts;
	rbpts.resize(degree);

	std::vector<double> alphas;
	alphas.resize(degree - 1);

	std::vector<double> errors;
	errors.resize(m,0.0);

	int nh = 0;

	updatedControlPoints.resize(degree);
	updatedControlPoints[0] = controlPoints[0];

	updatedKnotVector.resize(ph + 1);
	for (int i = 0; i <= ph; i++)
	{
		updatedKnotVector[i] = updatedKnotVector[0];
	}

	for (int i = 0; i <= static_cast<int>(degree); i++)
	{
		bpts[i] = controlPoints[i];
	}

	while (b < m)
	{
		int i = b;

		while (b < m && MathUtils::IsAlmostEqualTo(knotVector[b], knotVector[b + 1]))
		{
			b = b + 1;
		}
		mult = b - i + 1;
		mh += mult - 1;
		int oldr = r;
		r = degree - mult;

		int lbz = oldr > 0? (oldr + 2) / 2:1;

		if (r > 0)
		{
			double numer = knotVector[b] - knotVector[a];
			for (int k = degree; k >= mult; k--)
			{
				alphas[k - mult - 1] = numer / (knotVector[a + k] - knotVector[a]);
			}

			for (int j = 1; j <= r; j++)
			{
				int save = r - j;
				int s = mult + j;
				for (int k = degree; k >= s; k--)
				{
					bpts[k] = alphas[k - s] * bpts[k] + (1.0 - alphas[k - s]) * bpts[k - 1];
				}

				nextbpts[save] = bpts[degree];
			}
		}

		double maxError = ValidationUtils::ComputeMaxErrorOfBezierReduction(degree, bpts, rbpts);
		errors[a] += maxError;
		if (MathUtils::IsGreaterThan(errors[a], tol))
		{
			return false;
		}
		if (oldr > 0)
		{
			int first = kind;
			int last = kind;

			for (int k = 0; k < oldr; k++)
			{
				int i = first;
				int j = last;
				int kj = j - kind;

				while (j - i > k)
				{
					double alpha = (knotVector[a] - updatedKnotVector[i - 1]) / (knotVector[b] - updatedKnotVector[i - 1]);
					double beta = (knotVector[a] - updatedKnotVector[j - k - 1]) / (knotVector[b] - updatedKnotVector[j - k - 1]);
					updatedControlPoints[i - 1] = (updatedControlPoints[i - 1] - (1.0 - alpha) * updatedControlPoints[i - 2]) / alpha;
					rbpts[k] = (rbpts[kj] - beta * rbpts[kj + 1]) / (1.0 - beta);
					
					i = i + 1;
					j = j - 1;
					kj = kj - 1;
				}

				double Br;
				if (j - i < k)
				{
					Br = updatedControlPoints[i - 2].Distance(rbpts[kj + 1]);
				}
				else
				{
					double delta = (knotVector[a] - updatedKnotVector[i - 1]) / (knotVector[b] - updatedKnotVector[i - 1]);
					XYZW A = delta * rbpts[kj + 1] + (1.0 - delta) * updatedControlPoints[i - 2];
					Br = updatedControlPoints[i - 1].Distance(A);
				}
				
				int K = a + oldr - k;
				int q = (2 * degree - k + 1) / 2;
				int L = K - q;
				for (int ii = L; ii <= a; ii++)
				{
					errors[ii] += Br;
					if (MathUtils::IsGreaterThan(errors[ii], tol))
					{
						return false;
					}
				}

				first = first - 1;
				last = last - 1;
			}

			cind = i - 1;
		}

		if (a != degree)
		{
			for (int i = 0; i < ph - oldr; i++)
			{
				updatedKnotVector[kind] = knotVector[a];
				kind += 1;
			}
		}

		for (int i = lbz; i <= ph; i++)
		{
			updatedControlPoints[cind] = rbpts[i];
			cind += 1;
		}

		if (b < m)
		{
			for (int i = 0; i < r; i++)
			{
				bpts[i] = nextbpts[i];
			}
			for (int i = r; i <= static_cast<int>(degree); i++)
			{
				bpts[i] = controlPoints[b - degree + i];
			}

			a = b;
			b = b + 1;
		}
		else
		{
			for (int i = 0; i <= ph; i++)
			{
				updatedKnotVector[kind + i] = knotVector[b];
			}
		}
	}

	nh = mh - ph - 1;
	return true;
}

double LNLib::NurbsCurve::GetParamOnCurve(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, const XYZ& givenPoint)
{
	double minValue = Constants::MaxDistance;

	int maxIterations = 10;
	double paramT = Constants::DoubleEpsilon;
	double minParam = knotVector[0];
	double maxParam = knotVector[knotVector.size() - 1];

	int samples = static_cast<int>(controlPoints.size() * degree);
	double span = (maxParam - minParam) / (samples - 1);
	for (int i = 0; i < samples - 1; i++)
	{
		double currentU = minParam + span * i; 
		XYZ currentPoint;
		NurbsCurve::GetPointOnCurve(degree, knotVector, controlPoints, currentU, currentPoint);

		double nextU = minParam + span * (i + 1);
		XYZ nextPoint;
		NurbsCurve::GetPointOnCurve(degree, knotVector, controlPoints, nextU, nextPoint);

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
		std::vector<XYZ> derivatives;
		ComputeRationalCurveDerivatives(degree, knotVector, controlPoints, paramT, 2, derivatives);
		XYZ difference = derivatives[0] - givenPoint;
		double f = derivatives[1].DotProduct(difference);

		double condition1 = difference.Length();
		double condition2 = std::abs(f/ (derivatives[1].Length() * condition1));
		
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

void LNLib::NurbsCurve::CreateTransform(const std::vector<XYZW>& controlPoints, const Matrix4d& matrix, std::vector<XYZW>& transformedControlPoints)
{
	int size = static_cast<int>(controlPoints.size());
	transformedControlPoints.resize(size);

	Matrix4d tempMatrix = matrix;
	for (int i = 0; i < size; i++)
	{
		XYZW temp = controlPoints[i];
		transformedControlPoints[i] = XYZW(tempMatrix.OfPoint(temp.ToXYZ(true)),temp.GetW());
	}
}

void LNLib::NurbsCurve::ReverseKnotVector(const std::vector<double>& knotVector, std::vector<double> reversedKnotVector)
{
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
	int size = static_cast<int>(controlPoints.size());
	reversedControlPoints.resize(size);
	for (int i = 0; i < size; i++)
	{
		reversedControlPoints[i] = controlPoints[i];
	}
	std::reverse(reversedControlPoints.begin(), reversedControlPoints.end());
}

void LNLib::NurbsCurve::Reverse(const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& reversedKnotVector, std::vector<XYZW>& reversedControlPoints)
{
	ReverseKnotVector(knotVector, reversedKnotVector);
	ReverseControlPoints(controlPoints, reversedControlPoints);
}

bool LNLib::NurbsCurve::CreateArc(const XYZ& center, const XYZ& xAxis, const XYZ& yAxis, double startRad, double endRad, double xRadius, double yRadius, int& degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	XYZ nXTemp = xAxis;
	XYZ nX = nXTemp.Normalize();
	XYZ nYTemp = yAxis;
	XYZ nY = nYTemp.Normalize();

	if (endRad < startRad)
	{
		endRad = 2 * Constants::Pi + startRad;
	}
	double theta = endRad - startRad;
	
	int narcs = 0;
	if (MathUtils::IsLessThanOrEqual(theta, Constants::Pi / 2))
	{
		narcs = 1;
	}
	else
	{
		if (MathUtils::IsLessThanOrEqual(theta, Constants::Pi))
		{
			narcs = 2;
		}
		else if (MathUtils::IsLessThanOrEqual(theta, 3 * Constants::Pi / 2))
		{
			narcs = 3;
		}
		else
		{
			narcs = 4;
		}
	}
	double dtheta = theta / narcs;
	int n = 2 * narcs;

	knotVector.resize(n + 3);
	controlPoints.resize(n);

	double w1 = cos(dtheta / 2.0);
	XYZ P0 = center + xRadius * cos(startRad) * nX + yRadius * sin(startRad) * nY;
	XYZ T0 = -sin(startRad) * nX + cos(startRad) * nY;

	controlPoints[0] = XYZW(P0,1);
	int index = 0;
	double angle = startRad;
	for (int i = 1; i < narcs; i++)
	{
		angle = angle + dtheta;
		XYZ P2 = center + xRadius * cos(angle) * nX + yRadius * sin(angle) * nY;
		controlPoints[index + 2] = XYZW(P2, 1);
		XYZ T2 = -sin(angle) * nX + cos(angle) * nY;
		XYZ P1;
		double param0, param2 = 0.0;
		CurveCurveIntersectionType type =  Intersection::ComputeRays(P0, T0, P2, T2, param0, param2, P1);
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

	switch (narcs) {
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
	degree = 2;
	return true;
}

bool LNLib::NurbsCurve::CreateOneConicArc(const XYZ& start, const XYZ& startTangent, const XYZ& end, const XYZ& endTangent, const XYZ& pointOnConic, XYZ& projectPoint, double& projectPointWeight)
{
	XYZ sTemp = start;
	XYZ eTemp = end;

	double param0, param1 = 0.0;
	XYZ point = XYZ(0,0,0);
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
	else if(type == CurveCurveIntersectionType::Parallel)
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
	insertPointAtStartSide = (start + projectPointWeight * projectPoint) / (1 + projectPointWeight);
	insertPointAtEndSide = (projectPointWeight * projectPoint + end) / (1 + projectPointWeight);
	splitPoint = (insertPointAtStartSide + insertPointAtEndSide) / 2.0;
	insertWeight = sqrt((1 + projectPointWeight) / 2);
}

bool LNLib::NurbsCurve::CreateOpenConic(const XYZ& start, const XYZ& startTangent, const XYZ& end, const XYZ& endTangent, const XYZ& pointOnConic, int& degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
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

	controlPoints.resize(n + 1);
	knotVector.resize(j + 3);
	degree = 2;

	for (int i = 0; i < 3; i++)
	{
		knotVector[i] = 0.0;
		knotVector[i + j] = 1.0;
	}

	controlPoints[0] = XYZW(start,1.0);
	controlPoints[n] = XYZW(end, 1.0);

	if (nsegs == 1)
	{
		controlPoints[1] = XYZW(P1, w1);
		return true;
	}

	XYZ Q1,R1,S;
	double wqr = 0.0;
	SplitArc(start,P1,w1,end,Q1,S,R1,wqr);

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
		SplitArc(start,Q1,w1,S,HQ1,HS,HR1,wqr);
		controlPoints[2] = XYZW(HS, 1.0);
		controlPoints[1] = XYZW(HQ1, wqr);
		controlPoints[3] = XYZW(HR1, wqr);

		SplitArc(S,R1,w1,end,HQ1,HS,HR1,wqr);
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

void LNLib::NurbsCurve::GlobalInterpolation(unsigned int degree, const std::vector<XYZ>& throughPoints, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	int size = static_cast<int>(throughPoints.size());
	int n = size - 1;

	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	Interpolation::ComputeKnotVector(degree, static_cast<int>(throughPoints.size()), uk, knotVector);

	std::vector<std::vector<double>> A = Interpolation::MakeInterpolationMatrix(degree, size, uk, knotVector);
	std::vector<XYZ> tempControlPoints = Interpolation::GetSolvedMatrix(A, throughPoints);

	for (int i = 0; i < size; i++)
	{
		controlPoints[i] = XYZW(tempControlPoints[i], 1);
	}
}

void LNLib::NurbsCurve::GlobalInterpolationByTangents(unsigned int degree, const std::vector<XYZ>& throughPoints, const std::vector<XYZ>& tangents, const std::vector<int>& tangentsIndices, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	int pointsCount = static_cast<int>(throughPoints.size());
	int tangentsCount = static_cast<int>(tangents.size());
	int size = pointsCount + tangentsCount;
	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	Interpolation::ComputerKnotVectorForTangents(degree, uk, tangentsIndices, knotVector);

	std::vector<XYZ> innerTangents;
	innerTangents.resize(tangentsCount);

	double d = Interpolation::GetTotalChordLength(throughPoints);
	for (int i = 0; i < tangentsCount; i++)
	{
		XYZ t = tangents[i];
		XYZ tTemp = t.Normalize() * d;
		innerTangents[i] = tTemp;
	}

	std::vector<std::vector<double>> A;
	A.resize(size);
	for (int i = 0; i < size; i++)
	{
		A[i].resize(size);
	}

	std::vector<XYZ> b(size);

	int startRow = 1;
	int tangentStart = 0;

	if (tangentsIndices[0] == 0)
	{
		A[1][0] = -1;
		A[1][1] = 1;

		double y = (knotVector[degree + 1] - knotVector[0]) / degree;
		b[1] = y * innerTangents[0];
		startRow = 2;
		tangentStart = 1;
	}
	if (tangentsIndices[tangentsCount - 1] == pointsCount - 1)
	{
		A[size - 2][size - 2]     = -1;
		A[size - 2][size - 2 + 1] =  1;

		double y = (knotVector[knotVector.size() - 1] - knotVector[knotVector.size() - (degree + 2)]) / degree;
		b[b.size() - 2] = y * innerTangents[tangentsCount - 1];
	}

	int row = startRow;
	int tangentIndex = tangentStart;

	for (int i = 1; i < static_cast<int>(uk.size()) - 1; ++i)
	{
		int np = static_cast<int>(knotVector.size() - degree) - 2;
		int spanIndex = Polynomials::GetKnotSpanIndex(np, degree, uk[i], knotVector);

		if (tangentsIndices[tangentIndex] == i)
		{
			std::vector<std::vector<double>> derivatives;
			Polynomials::BasisFunctionsDerivatives(spanIndex, degree,uk[i],1,knotVector, derivatives);
			
			for (int r = 0; r < 2; r++)
			{
				for (int c = 0; c < static_cast<int>(degree) + 1; c++)
				{
					A[row + r][spanIndex - degree + c] = derivatives[r][c];
				}
			}
			
			b[row++] = throughPoints[i];
			b[row++] = innerTangents[tangentIndex++];
		}
		else
		{
			std::vector<double> basis;
			Polynomials::BasisFunctions(spanIndex, degree, uk[i], knotVector, basis);

			row++;
			for (int k = 0; k < static_cast<int>(degree); k++)
			{
				A[row][spanIndex - degree + k] = basis[k];
			}
		}
	}
	b[0] = throughPoints[0];
	b[b.size() - 1] = throughPoints[pointsCount - 1];
	A[0][0] = 1;
	A[size - 1][size - 1] = 1;

	std::vector<XYZ> tempControlPoints = Interpolation::GetSolvedMatrix(A, b);
	for (int i = 0; i < size; i++)
	{
		controlPoints[i] = XYZW(tempControlPoints[i], 1);
	}
}

void LNLib::NurbsCurve::LocalInterpolation(const std::vector<XYZ>& throughPoints, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	unsigned int degree = 3;

	int size = static_cast<int>(throughPoints.size());
	int n = size - 1;
	std::vector<XYZ> tangents;
	bool hasTangents = Interpolation::ComputerTangent(throughPoints, tangents);
	if (!hasTangents) return;

	std::vector<double> uk;
	uk.resize(n);
	uk[0] = 0;

	std::vector<XYZW> tempControlPoints;
	for (int k = 0; k < size; k++)
	{
		XYZ qk = throughPoints[k];

		XYZ t0 = tangents[k];
		XYZ t3 = tangents[k + 1];
		XYZ p0 = qk;
		XYZ p3 = throughPoints[k + 1];

		double a = 16 - (t0 + t3).Length();
		double b = 12 * (p3 - p0).DotProduct(t0 + t3);
		double c = -36 * (p3 - p0).SqrLength();

		double alpha = (-b + sqrt(b*b-4*a*c)) / (2 * a);

		XYZ pk0 = qk;
		XYZ pk1 = p0 + (1/3) * alpha * t0;
		XYZ pk2 = p3 - (1/3) * alpha * t3;

		uk[k + 1] = uk[k] + 3 * (pk1 - pk0).Length();

		tempControlPoints.emplace_back(XYZW(pk1,1));
		tempControlPoints.emplace_back(XYZW(pk2,1));
	}

	int kvSize = 2 * (degree + 1) + 2 * (n - 1);
	knotVector.resize(kvSize);
	for (int i = 0; i <= degree; i++)
	{
		knotVector[i] = 0;
		knotVector[kvSize - 1 - i] = 1;
	}

	for (int i = 1; i < n; i=i+2)
	{
		knotVector[degree + i] = knotVector[degree + (i + 1)] = uk[i] / uk[n];
	}

	int tSize = static_cast<int>(tempControlPoints.size());
	controlPoints.resize(tSize + 2);
	controlPoints[0] = XYZW(throughPoints[0], 1);
	controlPoints[tSize + 2 - 1] = XYZW(throughPoints[n], 1);
	for (int i = 0; i < tSize; i++)
	{
		controlPoints[i + 1] = tempControlPoints[i];
	}
}

bool LNLib::NurbsCurve::LeastSquaresApproximation(unsigned int degree, const std::vector<XYZ>& throughPoints, int controlPointsCount, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	int size = static_cast<int>(throughPoints.size());
	int n = controlPointsCount - 1;
	if (!ValidationUtils::IsInRange(n, degree, size - 1)) return false;

	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	int m = n + degree + 1;
	knotVector.resize(m,1.0);

	double d = (m + 1) / (n - degree + 1);
	for (int j = 0; j <= degree; j++)
		knotVector[j] = 0;

	for (int j = 1; j < n - degree; j++) 
	{
		knotVector[degree + j] = 0.0;
		for (int k = j; k < j + degree; ++k) 
		{
			int i = (int)(k * d);
			double a = k * d - i;
			int i2 = (int)((k - 1) * d);
			knotVector[degree + j] += a * uk[i2] + (1 - a) * uk[i];
		}
		knotVector[degree + j] /= degree;
	}

	std::vector<XYZ> R(n), rk(m);
	R[0] = throughPoints[0];
	R[n - 1] = throughPoints[m - 1];
	std::vector<std::vector<double>> N;
	N.resize(m);
	for (int i = 0; i < m; i++)
	{
		N[i].resize(n);
	}
	N[0][0] = 1.0;
	N[m - 1][n - 1] = 1.0;

	for (int i = 0; i < m; i++) 
	{
		int np = static_cast<int>(knotVector.size() - degree) - 2;
		int spanIndex = Polynomials::GetKnotSpanIndex(np, degree, uk[i], knotVector);
		std::vector<double> basis;
		Polynomials::BasisFunctions(spanIndex, degree, uk[i], knotVector, basis);
		for (int j = 0; j <= degree; j++) 
		{ 
			N[i][spanIndex - degree + j] = basis[j];
		}
		rk[i] = throughPoints[i] - N[i][0] * throughPoints[0] - N[i][n - 1] * throughPoints[m - 1];

	}

	for (int i = 0; i < n; i++) 
	{
		R[i] = XYZ(0,0,0);
		for (int j = 0; j < m; j++)
			R[i] += N[j][i] * rk[j];
		if (R[i].X() * R[i].X() < Constants::DoubleEpsilon &&
			R[i].Y() * R[i].Y() < Constants::DoubleEpsilon &&
			R[i].Z() * R[i].Z() < Constants::DoubleEpsilon)
		{
			return false;
		}
	}

	if (n - 2 > 0) 
	{
		std::vector<XYZ> B(n - 2);
		for (int i = 0; i < n - 2; i++)
		{
			B[i] = R[i + 1];
		}

		std::vector<std::vector<double>> Ns;
		Ns.resize(m - 2);
		for (int i = 0; i < m - 2; i++)
		{
			Ns[i].resize(n - 2);
		}
		for (int i = 1; i <= m - 2; i++)
		{
			for (int j = 1; j <= n - 2; j++)
			{
				Ns[i - 1][j - 1] = N[i][j];
			}
		}
		std::vector<std::vector<double>> Nst;
		MathUtils::Transpose(Ns, Nst);

		std::vector<XYZ> tempControlPoints = Interpolation::GetSolvedMatrix(MathUtils::MatrixMultiply(Nst,Ns), B);
		for (int i = 0; i < size; i++)
		{
			controlPoints[i] = XYZW(tempControlPoints[i], 1);
		}
	}

	controlPoints[0] = XYZW(throughPoints[0], 1);
	controlPoints[n - 1] = XYZW(throughPoints[m - 1], 1);
	return true;
}

double LNLib::NurbsCurve::ComputerRemoveKnotErrorBound(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, int removalIndex)
{
	int ord = static_cast<int>(degree + 1);
	int r = removalIndex;
	double u = knotVector[r];
	int s = Polynomials::GetKnotMultiplicity(u, knotVector);
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
	std::vector<int> S(knotSize,0);
	std::vector<int> Nl(knotSize);
	std::vector<int> Nr(knotSize,params.size() - 1);

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
			S[i] = Polynomials::GetKnotMultiplicity(knotVector[i],knotVector);
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
		int spanIndex = Polynomials::GetKnotSpanIndex(np, degree, uk[i], knotVector);
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
				spanIndex = Polynomials::GetKnotSpanIndex(np, degree, uk[i], tempNewU);
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
		P[i] = XYZW(throughPoints[i],1);
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

		bool result = LeastSquaresApproximation(deg + 1, throughPoints, P.size(), tU, tP);
		if (result)
		{
			U = tU;
			P = tP;

			for (int i = 0; i < size; i++)
			{
				double param = GetParamOnCurve(deg + 1, U, P, throughPoints[i]);
				uk[i] = param;
				XYZ p;
				GetPointOnCurve(deg + 1, U, P, param, p);
				error[i] = p.Distance(throughPoints[i]);
			}
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
		bool result = LeastSquaresApproximation(degree, throughPoints, P.size(), tU, tP);
		if (result)
		{
			for (int i = 0; i < size; i++)
			{
				double param = GetParamOnCurve(degree, U, P, throughPoints[i]);
				uk[i] = param;
				XYZ p;
				GetPointOnCurve(degree, U, P, param, p);
				error[i] = p.Distance(throughPoints[i]);
			}

			RemoveKnotsByGivenBound(degree, tU, tP, uk, error, maxError, Uh, Ph);
			knotVector = Uh;
			controlPoints = Ph;
			return;
		}
		else
		{
			knotVector = U;
			controlPoints = P;
			return;
		}
	}
}

void LNLib::NurbsCurve::ToUnclampCurve(unsigned int degree, std::vector<double>& knotVector, std::vector<XYZW>& controlPoints)
{
	int n = static_cast<int>(controlPoints.size() - 1);
	for (int i = 0; i < degree - 2; i++)
	{
		knotVector[degree - i - 1] = knotVector[degree - i] - (knotVector[n - i]);
		int k = degree - 1;
		for (int j = i; j >= 0; j--)
		{
			double alfa = (knotVector[degree] - knotVector[k]) / (knotVector[degree + j + 1] - knotVector[k]);
			controlPoints[j] = (controlPoints[j] - alfa * controlPoints[j + 1]) / (1.0 - alfa);
			k = k - 1;
		}
	}
	knotVector[0] = knotVector[1] - (knotVector[n + degree + 2] - knotVector[n - degree + 1]);
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
