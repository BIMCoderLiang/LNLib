/*
 * Author:
 * 2023/07/04 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Interpolation.h"
#include "XYZ.h"
#include "MathUtils.h"
#include "Intersection.h"
#include "Polynomials.h"

#include <algorithm>
#include <cmath>

namespace LNLib
{
	XYZ Getqk(const std::vector<XYZ>& throughPoints, int index)
	{
		return throughPoints[index] - throughPoints[index - 1];
	}

	double Getak(const XYZ& qk_1, const XYZ& qk, const XYZ& qk1, const XYZ& qk2)
	{
		return (qk_1.CrossProduct(qk)).Length() / ((qk_1.CrossProduct(qk).Length()) + (qk1.CrossProduct(qk2)).Length());
	}

	XYZ GetTk(const XYZ& qk_1, const XYZ& qk, const XYZ& qk1, const XYZ& qk2)
	{
		double ak = Getak(qk_1, qk, qk1, qk2);
		return ((1 - ak) * qk + ak * qk1).Normalize();
	}
}

double LNLib::Interpolation::GetTotalChordLength(const std::vector<XYZ>& throughPoints)
{
	int n = throughPoints.size() - 1;
	double length = 0.0;
	for (int i = 1; i <= n; i++)
	{
		length += throughPoints[i].Distance(throughPoints[i - 1]);
	}
	return length;
}

std::vector<double> LNLib::Interpolation::GetChordParameterization(const std::vector<XYZ>& throughPoints)
{
	int size = throughPoints.size();
	int n = size - 1;

	std::vector<double> uk(size,0.0);
	uk[n] = 1.0;

	double d = GetTotalChordLength(throughPoints);
	for (int i = 1; i <= n - 1; i++)
	{
		uk[i] = uk[i - 1] + (throughPoints[i].Distance(throughPoints[i - 1])) / d;
	}
	return uk;
}

double LNLib::Interpolation::GetCentripetalLength(const std::vector<XYZ>& throughPoints)
{
	int size = throughPoints.size();
	int n = size - 1;

	double length = 0.0;
	for (int i = 1; i <= n; i++)
	{
		length += std::sqrt(throughPoints[i].Distance(throughPoints[i - 1]));
	}
	return length;
}

std::vector<double> LNLib::Interpolation::GetCentripetalParameterization(const std::vector<XYZ>& throughPoints)
{
	int size = throughPoints.size();
	int n = size - 1;

	std::vector<double> uk(size, 0.0);
	uk[n] = 1.0;

	double d = GetCentripetalLength(throughPoints);
	for (int i = 1; i <= n - 1; i++)
	{
		uk[i] = uk[i - 1] + std::sqrt(throughPoints[i].Distance(throughPoints[i - 1])) / d;
	}
	return uk;
}

std::vector<double> LNLib::Interpolation::AverageKnotVector(int degree, const std::vector<double>& params)
{
	
	std::vector<double> uk = params;
	int size = params.size();
	int n = size - 1;
	int m = n + degree + 1;

	std::vector<double> knotVector(m + 1, 0.0);
	for (int i = m - degree; i <= m; i++)
	{
		knotVector[i] = 1.0;
	}

	for (int j = 1; j <= n - degree; j++)
	{
		double sum = 0.0;
		for (int i = j; i <= j + degree - 1; i++)
		{
			sum += uk[i];
		}
		knotVector[j + degree] = (1.0 / degree) * sum;
	}
	return knotVector;
}

std::vector<double> LNLib::Interpolation::ComputeKnotVector(int degree, int controlPointsCount, const std::vector<double>& params)
{
	int m = params.size() - 1;
	int n = controlPointsCount - 1;
	int nn = n + degree + 2;

	std::vector<double>  knotVector(nn,0.0);
	double d = (double)(m + 1) / (double)(n - degree + 1);
	for (int j = 1; j <= n - degree; j++)
	{
		int i = floor(j * d);
		double alpha = (j * d) - i;
		knotVector[degree + j] = (1.0 - alpha) * params[i - 1] + (alpha * params[i]);
	}
	for (int i = 0; i < nn; i++)
	{
		if (i <= degree)
		{
			knotVector[i] = 0.0;
		}
		else if (i >= nn - 1 - degree)
		{
			knotVector[i] = 1.0;
		}
	}
	return knotVector;
}

bool LNLib::Interpolation::ComputerWeightForRationalQuadraticInterpolation(const XYZ& startPoint, const XYZ& middleControlPoint, const XYZ& endPoint, double& weight)
{
	XYZ SM = middleControlPoint - startPoint;
	XYZ EM = middleControlPoint - endPoint;
	XYZ SE = endPoint - startPoint;

	if (SM.Normalize().IsAlmostEqualTo(EM.Normalize()) ||
		SM.Normalize().IsAlmostEqualTo(EM.Normalize().Negative()))
	{
		weight = 1.0;
		return true;
	}

	if (MathUtils::IsAlmostEqualTo(SM.Length(), EM.Length()))
	{
		weight = SM.DotProduct(SE) / (SM.Length() * SE.Length());
		return true;
	}

	XYZ M = 0.5 * (startPoint + endPoint);
	XYZ MR = middleControlPoint - M;
	double ratio = SM.Length() / SE.Length();
	double frac = 1.0 / (ratio + 1.0);
	XYZ D = endPoint + frac * EM;
	XYZ SD = D - startPoint;
	XYZ S1;

	double alf1 = 0.0;
	double alf2 = 0.0;

	auto type = Intersection::ComputeRays(startPoint, SD.Normalize(), M, MR.Normalize(), alf1, alf2, S1);
	if (type != CurveCurveIntersectionType::Intersecting)
	{
		return false;
	}

	ratio = EM.Length() / SE.Length();
	frac = 1.0 / (ratio + 1.0);
	D = startPoint + frac * SM;
	XYZ ED = D - endPoint;
	XYZ S2;
	type = Intersection::ComputeRays(endPoint, ED.Normalize(), M, MR.Normalize(), alf1, alf2, S2);
	if (type != CurveCurveIntersectionType::Intersecting)
	{
		return false;
	}

	XYZ S = 0.5 * (S1 + S2);
	XYZ MS = S - M;
	double s = MS.Length() / MR.Length();
	weight = s / (1.0 - s);
	return true;
}

bool LNLib::Interpolation::GetSurfaceMeshParameterization(const std::vector<std::vector<XYZ>>& throughPoints, std::vector<double>& paramsU, std::vector<double>& paramsV)
{
	int n = throughPoints.size();
	int m = throughPoints[0].size();

	std::vector<double> cds(std::max(n, m),0.0);
	paramsU.resize(n,0.0);
	paramsV.resize(m,0.0);

	int num = m;
	
	for (int l = 0; l < m; l++)
	{
		double total = 0.0;
		for (int k = 1; k < n; k++)
		{
			cds[k] = throughPoints[k][l].Distance(throughPoints[k - 1][l]);
			total += cds[k];
		}

		if (MathUtils::IsAlmostEqualTo(total, 0.0))
		{
			num--;
		}
		else
		{
			double d = 0.0;
			for (int k = 1; k < n; k++)
			{
				d += cds[k];
				paramsU[k] = paramsU[k] + d / total;
			}
		}
	}
	if (num == 0)
	{
		return false;
	}

	for (int k = 1; k < n - 1; k++)
	{
		paramsU[k] = paramsU[k] / num;
	}
	paramsU[n - 1] = 1.0;

	num = n;

	for (int k = 0; k < n; k++) 
	{
		double total = 0.0;
		for (int l = 1; l < m; l++) 
		{
			cds[l] = throughPoints[k][l].Distance(throughPoints[k][l-1]);
			total += cds[l];
		}
		if (MathUtils::IsAlmostEqualTo(total, 0.0))
		{
			num--;
		}	
		else 
		{
			double d = 0.0;
			for (int l = 1; l < m; l++) 
			{
				d += cds[l];
				paramsV[l] += d / total;
			}
		}
	}

	if (num == 0)
	{
		return false;
	}
	for (int l = 1; l < m - 1; l++)
	{
		paramsV[l] = paramsV[l] / num;
	}
	paramsV[m - 1] = 1.0;

	return true;
}

bool LNLib::Interpolation::ComputeTangent(const std::vector<XYZ>& throughPoints, std::vector<XYZ>& tangents)
{
	int size = throughPoints.size();
	if (size < 5) return false;

	tangents.resize(size);
	for (int k = 2; k < size - 2 ; k++)
	{
		LNLib::XYZ qk_1 = Getqk(throughPoints, k - 1);
		LNLib::XYZ qk = Getqk(throughPoints, k);
		LNLib::XYZ qk1 = Getqk(throughPoints, k+1);
		LNLib::XYZ qk2 = Getqk(throughPoints, k+2);

		tangents[k] = GetTk(qk_1, qk, qk1, qk2);
	}

	int n = size - 1;
	LNLib::XYZ q0 = 2* Getqk(throughPoints, 1) - Getqk(throughPoints, 2);
	LNLib::XYZ q_1 = 2 * q0 - Getqk(throughPoints, 1);
	LNLib::XYZ qn1 = 2 * Getqk(throughPoints, n) - Getqk(throughPoints, n-1);
	LNLib::XYZ qn2 = 2 * qn1 - Getqk(throughPoints, n);

	tangents[0] = GetTk(q_1, q0, Getqk(throughPoints, 1), Getqk(throughPoints, 2));
	tangents[1] = GetTk(q0, Getqk(throughPoints, 1), Getqk(throughPoints, 2), Getqk(throughPoints, 3));
	tangents[n-1] = GetTk(Getqk(throughPoints, n - 2), Getqk(throughPoints, n-1), Getqk(throughPoints, n),qn1);
	tangents[n] = GetTk(Getqk(throughPoints, n-1), Getqk(throughPoints, n),qn1,qn2);

	return true;
}

std::vector<LNLib::XYZ> LNLib::Interpolation::ComputeTangent(const std::vector<XYZ>& throughPoints)
{
	int size = throughPoints.size();
	std::vector<XYZ> tangents(size, XYZ(0, 0, 0));
	std::vector<XYZ> qq(size, XYZ(0, 0, 0));
	std::vector<double> delta(size, 0.0);

	auto params = GetChordParameterization(throughPoints);
	for (int i = 1; i < size; i++)
	{
		delta[i] = params[i] - params[i - 1];
		qq[i] = throughPoints[i] - throughPoints[i - 1];
	}
	for (int i = 1; i < size - 1; i++)
	{
		double a = delta[i] / (delta[i] + delta[i + 1]);
		tangents[i] = ((1 - a) * qq[i] + a * qq[i + 1]).Normalize();
	}

	tangents[0] = (2 * qq[1] / delta[1] - tangents[1]).Normalize();
	tangents[tangents.size() - 1] = (2 * qq[qq.size() - 1] / delta[delta.size() - 1] - tangents[tangents.size() - 2]).Normalize();
	return tangents;
}
