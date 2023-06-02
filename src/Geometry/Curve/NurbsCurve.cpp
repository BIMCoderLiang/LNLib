#include "NurbsCurve.h"
#include "Constants.h"
#include "Polynomials.h"
#include "XYZ.h"
#include "XYZW.h"
#include "MathUtils.h"
#include "BsplineCurve.h"
#include "ValidationUtils.h"
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
	int np = static_cast<int>(knotVector.size() - degree) - 2;
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
	int n = static_cast<int>(knotVector.size() - degree) - 2;

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
	int n = static_cast<int>(knotVector.size() - degree - 2);
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
	int n = static_cast<int>(knotVector.size() - degree - 2);
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

void LNLib::NurbsCurve::RemoveKnot(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double removeKnot, unsigned int removeIndex, unsigned int times, std::vector<double>& restKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);

	int n = static_cast<int>(knotVector.size() - degree - 2);
	int m = n + degree + 1;

	int order = degree + 1;
	int originMultiplicity = Polynomials::GetKnotMultiplicity(removeKnot, knotVector);
	double fout = (2 * removeIndex - originMultiplicity - degree) / 2;
	int last = removeIndex - originMultiplicity;
	int first = removeIndex - degree;

	for (int i = 0; i < static_cast<int>(controlPoints.size()); i++)
	{
		updatedControlPoints[i] = controlPoints[i];
	}

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
		int remflag = 0;
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
				remflag = 1;
			}
		}
		else
		{
			double alphai = (removeKnot - knotVector[i]) / (knotVector[i + order + t] - knotVector[i]);
			if (MathUtils::IsLessThanOrEqual(controlPoints[i].Distance(alphai * temp[ii + t + 1] + (1.0 * alphai) * temp[ii - 1]), tol))
			{
				remflag = 1;
			}
		}

		if (remflag == 0)
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
	for (int k = removeIndex + 1; k <= m; k++)
	{
		restKnotVector[k - t] = knotVector[k];
	}

	int j = static_cast<int>(fout);
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
}

void LNLib::NurbsCurve::ElevateDegree(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, unsigned int times,  std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	int n = static_cast<int>(knotVector.size() - degree - 2);
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
			bezalfs[i][j] = inv * MathUtils::Binomial(times, i - j);
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
	int b = degree + 1;
	int cind = 1;
	double ua = knotVector[0];
	updatedControlPoints[0] = controlPoints[0];

	for (int i = 0; i <= ph; i++)
	{
		updatedKnotVector[i] = ua;
	}

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
		r = degree - mul;

		int lbz = 0;

		if (oldr > 0)
		{
			lbz = (oldr + 2) / 2;
		}
		else
		{
			lbz = 1;
		}

		double rbz = 0.0;
		if (r > 0)
		{
			rbz = ph - (r + 1) / 2;
		}
		else
		{
			rbz = ph;
		}

		if (r > 0)
		{
			double numer = ub - ua;
			std::vector<double> alfs;
			alfs.resize(degree - 1);
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

		nh = mh - ph - 1;
	}
}

int LNLib::NurbsCurve::ReduceDegree(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, std::vector<double>& updatedKnotVector, std::vector<XYZW>& updatedControlPoints)
{
	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);

	int ph = degree - 1;
	int mh = ph;

	int kind = ph + 1;
	int r = -1;
	int a = degree;

	int b = degree + 1;
	int cind = 1;
	int mult = degree;

	int n = static_cast<int>(knotVector.size() - degree - 2);
	int m = n + degree + 1;

	std::vector<XYZW> bpts;
	bpts.resize(degree + 1);

	std::vector<XYZW> nextbpts;
	nextbpts.resize(degree - 1);

	std::vector<XYZW> rbpts;
	rbpts.resize(degree);

	std::vector<double> alphas;
	alphas.resize(degree - 1);

	std::vector<double> error;
	error.resize(m);

	int nh = 0;

	updatedControlPoints.resize(degree);
	updatedControlPoints[0] = controlPoints[0];

	for (int i = 0; i <= ph; i++)
	{
		updatedKnotVector[i] = updatedKnotVector[0];
	}

	for (int i = 0; i <= static_cast<int>(degree); i++)
	{
		bpts[i] = controlPoints[i];
	}

	for (int i = 0; i < m; i++)
	{
		error[i] = 0.0;
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

		int lbz = 0;
		if (oldr > 0)
		{
			int lbz = (oldr + 2) / 2;
		}
		else
		{
			lbz = 1;
		}

		if (r > 0)
		{
			double numer = knotVector[b] - knotVector[a];
			for (int k = degree; k >= mult; k--)
			{
				alphas[k - mult - 1] = numer / (knotVector[a + k] - knotVector[a]);
			}

			for (int j = 1; j < r; j++)
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
		error[a] += maxError;
		if (MathUtils::IsGreaterThan(error[a], tol))
		{
			return 1;
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
					error[ii] += Br;
					if (MathUtils::IsGreaterThan(error[ii], tol))
					{
						return 1;
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
	return 0;
}
