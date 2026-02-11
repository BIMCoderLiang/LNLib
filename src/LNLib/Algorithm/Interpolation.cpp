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
	size_t size = throughPoints.size();

	if (size == 0) {
		return {};
	}
	if (size == 1) {
		return { 0.0 };
	}
	if (size == 2) {
		return { 0.0, 1.0 };
	}

	std::vector<double> segmentLengths;
	segmentLengths.reserve(size - 1);
	double totalLength = 0.0;

	for (size_t i = 1; i < size; ++i) {
		double segLen = throughPoints[i].Distance(throughPoints[i - 1]);
		segmentLengths.push_back(segLen);
		totalLength += segLen;
	}

	if (MathUtils::IsAlmostEqualTo(totalLength, 0.0)) {
		std::vector<double> uk(size);
		for (size_t i = 0; i < size; ++i) {
			uk[i] = static_cast<double>(i) / static_cast<double>(size - 1);
		}
		return uk;
	}

	std::vector<double> uk(size);
	uk[0] = 0.0;
	uk[size - 1] = 1.0;

	double accumulated = 0.0;
	for (size_t i = 1; i < size - 1; ++i) {
		accumulated += segmentLengths[i - 1];
		uk[i] = accumulated / totalLength;
	}

	return uk;
}

double LNLib::Interpolation::GetCentripetalLength(const std::vector<XYZ>& throughPoints)
{
	if (throughPoints.size() < 2)
		return 0.0;

	double length = 0.0;
	for (size_t i = 1; i < throughPoints.size(); ++i) {
		double dist = throughPoints[i].Distance(throughPoints[i - 1]);
		length += std::sqrt(dist);
	}
	return length;
}

std::vector<double> LNLib::Interpolation::GetCentripetalParameterization(const std::vector<XYZ>& throughPoints)
{
	const size_t size = throughPoints.size();

	if (size == 0) {
		return {};
	}
	if (size == 1) {
		return { 0.0 };
	}
	if (size == 2) {
		return { 0.0, 1.0 };
	}

	std::vector<double> segmentLengths;
	segmentLengths.reserve(size - 1);
	double totalLength = 0.0;

	for (size_t i = 1; i < size; ++i) {
		double euclideanDist = throughPoints[i].Distance(throughPoints[i - 1]);
		double centripetalSeg = std::sqrt(euclideanDist);
		segmentLengths.push_back(centripetalSeg);
		totalLength += centripetalSeg;
	}

	if (MathUtils::IsAlmostEqualTo(totalLength, 0.0)) {
		std::vector<double> uk(size);
		for (size_t i = 0; i < size; ++i) {
			uk[i] = static_cast<double>(i) / static_cast<double>(size - 1);
		}
		return uk;
	}

	std::vector<double> uk(size);
	uk[0] = 0.0;
	uk[size - 1] = 1.0;

	double accumulated = 0.0;
	for (size_t i = 1; i < size - 1; ++i) {
		accumulated += segmentLengths[i - 1];
		uk[i] = accumulated / totalLength;
	}

	return uk;
}

std::vector<double> LNLib::Interpolation::AverageKnotVector(int degree, const std::vector<double>& params)
{
	int n = static_cast<int>(params.size()) - 1;
	int m = n + degree + 1;
	std::vector<double> knotVector(m + 1, 0.0);

	for (int i = m - degree; i <= m; ++i) {
		knotVector[i] = 1.0;
	}

	if (n <= degree) {
		return knotVector;
	}

	double windowSum = 0.0;
	for (int i = 1; i <= degree; ++i) {
		windowSum += params[i];
	}
	knotVector[degree + 1] = windowSum / static_cast<double>(degree);

	for (int j = 2; j <= n - degree; ++j) {
		windowSum -= params[j - 1];
		windowSum += params[j + degree - 1];
		knotVector[degree + j] = windowSum / static_cast<double>(degree);
	}

	return knotVector;
}

std::vector<double> LNLib::Interpolation::ComputeKnotVector(int degree, int controlPointsCount, const std::vector<double>& params)
{
	int m = params.size() - 1;
	int n = controlPointsCount - 1;
	int nn = n + degree + 2;

	std::vector<double>  knotVector(nn, 0.0);
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

	std::vector<double> cds(std::max(n, m), 0.0);
	paramsU.resize(n, 0.0);
	paramsV.resize(m, 0.0);

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
			cds[l] = throughPoints[k][l].Distance(throughPoints[k][l - 1]);
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
	for (int k = 2; k < size - 2; k++)
	{
		LNLib::XYZ qk_1 = Getqk(throughPoints, k - 1);
		LNLib::XYZ qk = Getqk(throughPoints, k);
		LNLib::XYZ qk1 = Getqk(throughPoints, k + 1);
		LNLib::XYZ qk2 = Getqk(throughPoints, k + 2);

		tangents[k] = GetTk(qk_1, qk, qk1, qk2);
	}

	int n = size - 1;
	LNLib::XYZ q0 = 2 * Getqk(throughPoints, 1) - Getqk(throughPoints, 2);
	LNLib::XYZ q_1 = 2 * q0 - Getqk(throughPoints, 1);
	LNLib::XYZ qn1 = 2 * Getqk(throughPoints, n) - Getqk(throughPoints, n - 1);
	LNLib::XYZ qn2 = 2 * qn1 - Getqk(throughPoints, n);

	tangents[0] = GetTk(q_1, q0, Getqk(throughPoints, 1), Getqk(throughPoints, 2));
	tangents[1] = GetTk(q0, Getqk(throughPoints, 1), Getqk(throughPoints, 2), Getqk(throughPoints, 3));
	tangents[n - 1] = GetTk(Getqk(throughPoints, n - 2), Getqk(throughPoints, n - 1), Getqk(throughPoints, n), qn1);
	tangents[n] = GetTk(Getqk(throughPoints, n - 1), Getqk(throughPoints, n), qn1, qn2);

	return true;
}

std::vector<LNLib::XYZ> LNLib::Interpolation::ComputeTangent(const std::vector<XYZ>& throughPoints)
{
	size_t size = throughPoints.size();

	if (size == 0) return {};
	if (size == 1) return { XYZ(0, 0, 0) };
	if (size == 2) {
		XYZ dir = throughPoints[1] - throughPoints[0];
		return { dir, dir };
	}

	auto params = GetChordParameterization(throughPoints);

	std::vector<double> du(size - 1);
	std::vector<XYZ> dP(size - 1);
	for (size_t i = 0; i < size - 1; ++i) {
		du[i] = params[i + 1] - params[i];
		dP[i] = throughPoints[i + 1] - throughPoints[i];
	}

	std::vector<XYZ> tangents(size);

	for (size_t i = 1; i < size - 1; ++i) {
		double w1 = du[i];
		double w2 = du[i - 1];
		if (MathUtils::IsAlmostEqualTo(w1, 0.0) &&
			MathUtils::IsAlmostEqualTo(w2, 0.0)) {
			tangents[i] = XYZ(0, 0, 0);
		}
		else if (MathUtils::IsAlmostEqualTo(w1, 0.0)) {
			tangents[i] = dP[i] / du[i];
		}
		else if (MathUtils::IsAlmostEqualTo(w2, 0.0)) {
			tangents[i] = dP[i - 1] / du[i - 1];
		}
		else {
			double a = w2 / (w1 + w2);
			tangents[i] = a * (dP[i - 1] / du[i - 1]) + (1 - a) * (dP[i] / du[i]);
		}
	}

	if (MathUtils::IsGreaterThan(du[0], 0.0)) {
		XYZ deriv1 = dP[0] / du[0];
		tangents[0] = 2.0 * deriv1 - tangents[1];
	}
	else {
		tangents[0] = tangents[1];
	}

	size_t last = size - 1;
	if (MathUtils::IsGreaterThan(du[last - 1], 0.0)) {
		XYZ derivN = dP[last - 1] / du[last - 1];
		tangents[last] = 2.0 * derivN - tangents[last - 1];
	}
	else {
		tangents[last] = tangents[last - 1];
	}

	return tangents;
}
