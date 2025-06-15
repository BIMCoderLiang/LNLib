/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
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
#include "KnotVectorUtils.h"
#include "ControlPointsUtils.h"
#include "Interpolation.h"
#include "Integrator.h"
#include "LNLibExceptions.h"
#include "LNObject.h"
#include <vector>
#include <set>
#include <random>
#include <algorithm>
#include <cmath>

namespace LNLib
{
	class FirstDerivativeLengthFunction : public IntegrationFunction
	{
		double operator()(double parameter, void* customData)
		{
			LN_NurbsCurve* curve = (LN_NurbsCurve*)customData;
			return NurbsCurve::ComputeRationalCurveDerivatives(*curve, 1, parameter)[1].Length();
		}
	};

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

	double CalculateLengthBySimpson(FirstDerivativeLengthFunction function, const LN_NurbsCurve& curve, double start, double end, double simpson, double tolearance)
	{
		double length = 0.0;
		double m = (start + end) / 2.0;
		double left = Integrator::Simpson(function, (void*)& curve, start, m);
		double right = Integrator::Simpson(function, (void*)& curve, m, end);

		double differ = left + right - simpson;
		if (MathUtils::IsAlmostEqualTo(differ, 0.0) || 
			MathUtils::IsLessThan(abs(differ) / 10.0, tolearance))
		{
			length = left + right + differ / 10.0;
		}
		else
		{
			length = CalculateLengthBySimpson(function, curve, start, m, left, tolearance / 2.0) + CalculateLengthBySimpson(function, curve, m, end, right, tolearance / 2.0);
		}
		return length;
	}

	double GetParamByLength(const LN_NurbsCurve& curve, double start, double end, double givenLength, IntegratorType type)
	{
		double middle = (start + end) / 2.0;

		LN_NurbsCurve left;
		LN_NurbsCurve right;
		bool isSplited = NurbsCurve::SplitAt(curve, middle, left, right);
		if (!isSplited) return start;
		double length = NurbsCurve::ApproximateLength(left, type);
		if (MathUtils::IsAlmostEqualTo(length, givenLength, Constants::DistanceEpsilon))
		{
			return middle;
		}
		else if (MathUtils::IsGreaterThan(length, givenLength, Constants::DistanceEpsilon))
		{
			end = middle;
			return GetParamByLength(curve, start, end, givenLength, type);
		}
		else if (MathUtils::IsLessThan(length, givenLength, Constants::DistanceEpsilon))
		{
			start = middle;
			return GetParamByLength(curve, start, end, givenLength, type);
		}
		return middle;
	}

	void TessellateCore(const LN_NurbsCurve& curve, double start, double end, std::vector<double>& parameters)
	{
		double half = 0.5;
		double mid = start + (end - start) * half;
		
		XYZ startPoint = NurbsCurve::GetPointOnCurve(curve, start);
		XYZ endPoint = NurbsCurve::GetPointOnCurve(curve, end);

		XYZ midPoint = NurbsCurve::GetPointOnCurve(curve, mid);
		XYZ segMidPoint = startPoint + (endPoint - startPoint) * half;

		if (MathUtils::IsLessThanOrEqual(midPoint.Distance(endPoint), Constants::DistanceEpsilon))
		{
			parameters.emplace_back(mid);
		}
		else
		{
			double distance = midPoint.Distance(segMidPoint);
			bool condition1 = MathUtils::IsLessThanOrEqual(distance, Constants::DistanceEpsilon);

			XYZ startTangent = NurbsCurve::ComputeRationalCurveDerivatives(curve, 1, start)[1];
			XYZ endTangent = NurbsCurve::ComputeRationalCurveDerivatives(curve, 1, end)[1];
			double angle = startTangent.AngleTo(endTangent);
			bool condition2 = MathUtils::IsLessThanOrEqual(angle, Constants::AngleEpsilon);

			bool isIgnoreCurvature = false;
			bool condition3 = false;
			double halfChordLength = startPoint.Distance(endPoint) / 2.0;
			double curvature = NurbsCurve::Curvature(curve, mid);
			if (MathUtils::IsAlmostEqualTo(curvature, 0.0))
			{
				isIgnoreCurvature = true;
			}
			else
			{
				double radius = 1.0 / NurbsCurve::Curvature(curve, mid);
				double chordHeight = radius - std::sqrt(radius * radius - halfChordLength * halfChordLength);
				condition3 = MathUtils::IsLessThanOrEqual(chordHeight, Constants::DistanceEpsilon);
			}
			bool condition = isIgnoreCurvature ? (condition1 && condition2) : (condition1 && condition2 && condition3);
			if (condition)
			{
				parameters.emplace_back(mid);
			}
			else
			{
				TessellateCore(curve, start, mid, parameters);
				TessellateCore(curve, mid, end, parameters);
			}
		}
	}

	void OffsetByTillerAndHanson(const LN_NurbsCurve& curve, double offset, LN_NurbsCurve& result)
	{
		int degree = curve.Degree;
		std::vector<double> knotVector = curve.KnotVector;
		std::vector<XYZW> controlPoints = curve.ControlPoints;

		if (MathUtils::IsAlmostEqualTo(offset, 0.0))
		{
			result = curve;
			return;
		}

		std::vector<XYZ> tempPoints;

		for (int i = 0; i < controlPoints.size() - 1; i++)
		{
			XYZ currentPoint0 = controlPoints[i].ToXYZ(true);
			double param0 = NurbsCurve::GetParamOnCurve(curve, currentPoint0);

			XYZ currentPoint1 = controlPoints[i + 1].ToXYZ(true);
			double param1 = NurbsCurve::GetParamOnCurve(curve, currentPoint1);

			XYZ newPoint0 = currentPoint0 + offset * NurbsCurve::Normal(curve, CurveNormal::Normal, param0);
			XYZ newPoint1 = currentPoint1 + offset * NurbsCurve::Normal(curve, CurveNormal::Normal, param1);

			tempPoints.emplace_back(newPoint0);
			tempPoints.emplace_back(newPoint1);
		}

		std::vector<XYZW> newControlPoints;
		newControlPoints.reserve(controlPoints.size());
		newControlPoints.emplace_back(XYZW(tempPoints[0], 1));

		for (int i = 0; i < tempPoints.size() - 2; i = i + 2)
		{
			XYZ currentP0 = tempPoints[i];
			XYZ currentP1 = tempPoints[i + 1];

			XYZ nextP0 = tempPoints[i + 2];
			XYZ nextP1 = tempPoints[i + 3];

			double param0, param1;
			XYZ result;
			CurveCurveIntersectionType type = Intersection::ComputeRays(currentP0, currentP1 - currentP0, nextP1, nextP0 - nextP1, param0, param1, result);
			if (type != CurveCurveIntersectionType::Intersecting)
			{
				return;
			}
			newControlPoints.emplace_back(XYZW(result, 1));
		}

		newControlPoints.emplace_back(XYZW(tempPoints[tempPoints.size() - 1], 1));

		result = curve;
		result.ControlPoints = newControlPoints;
	}

	void OffsetByPieglAndTiller(const LN_NurbsCurve& curve, double offset, LN_NurbsCurve& result)
	{
		int degree = curve.Degree;
		std::vector<double> knotVector = curve.KnotVector;
		std::vector<XYZW> controlPoints = curve.ControlPoints;

		if (MathUtils::IsAlmostEqualTo(offset, 0.0))
		{
			result = curve;
			return;
		}

		double minParam = knotVector[0];
		double maxParam = knotVector[knotVector.size() - 1];

		int count = 2 * (degree + 1);
		int n = 0;
		for (int i = 0; i < knotVector.size(); i++)
		{
			std::vector<XYZ> ders = NurbsCurve::ComputeRationalCurveDerivatives(curve, 2, knotVector[i]);
			if (ders.size() > 2)
			{
				XYZ d2 = NurbsCurve::ComputeRationalCurveDerivatives(curve, 2, knotVector[i])[2];
				int l = (int)(sqrt(1.0 / Constants::DistanceEpsilon) * sqrt(d2.Length() / 8.0));
				if (l > n)
				{
					n = l;
				}
			}
		}
		if (n > count)
		{
			count = n;
		}
		double step = (maxParam - minParam) / (double)(count - 1);
		std::vector<XYZ> samplePoints;
		samplePoints.reserve(count);
		for (int j = 0; j < count; j++)
		{
			double bparam = minParam + step * j;
			XYZ point = NurbsCurve::GetPointOnCurve(curve, bparam);
			XYZ newPoint = point + offset * NurbsCurve::Normal(curve, CurveNormal::Normal, bparam);
			samplePoints.emplace_back(newPoint);
		}

		LN_NurbsCurve temp;
		NurbsCurve::GlobalInterpolation(degree, samplePoints, temp);
		NurbsCurve::RemoveExcessiveKnots(temp, result);
	}
}

void LNLib::NurbsCurve::Check(const LN_NurbsCurve& curve)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must be greater than zero.");
	VALIDATE_ARGUMENT(knotVector.size() > 0, "knotVector", "KnotVector size must be greater than zero.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidKnotVector(knotVector), "knotVector", "KnotVector must be a nondecreasing sequence of real numbers.");
	VALIDATE_ARGUMENT(controlPoints.size() > 0, "controlPoints", "ControlPoints must contain one point at least.");
	VALIDATE_ARGUMENT(ValidationUtils::IsValidNurbs(degree, knotVector.size(), controlPoints.size()), "controlPoints", "Arguments must be fit: m = n + p + 1");
}


LNLib::XYZ LNLib::NurbsCurve::GetPointOnCurve(const LN_NurbsCurve& curve, double paramT)
{
	XYZW weightPoint = BsplineCurve::GetPointOnCurve(curve, paramT);
	return weightPoint.ToXYZ(true);
}

std::vector<LNLib::XYZ> LNLib::NurbsCurve::ComputeRationalCurveDerivatives(const LN_NurbsCurve& curve, int derivative, double paramT)
{
	std::vector<LNLib::XYZ> derivatives(derivative + 1);
	std::vector<XYZW> ders = BsplineCurve::ComputeDerivatives(curve, derivative, paramT);

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

bool LNLib::NurbsCurve::CanComputerDerivative(const LN_NurbsCurve& curve, double paramT)
{
	// Treat boundary parameters as derivable.
	if(MathUtils::IsAlmostEqualTo(paramT, curve.KnotVector[0]) ||
	   MathUtils::IsAlmostEqualTo(paramT, curve.KnotVector.back()))
	{
		return true;
	}

	auto pt = GetPointOnCurve(curve, paramT);
	XYZ left = (GetPointOnCurve(curve, paramT - Constants::DoubleEpsilon) - pt) / -Constants::DoubleEpsilon;
	XYZ right = (GetPointOnCurve(curve, paramT + Constants::DoubleEpsilon) - pt) / Constants::DoubleEpsilon;

	return left.IsAlmostEqualTo(right);
}

double LNLib::NurbsCurve::Curvature(const LN_NurbsCurve& curve, double paramT)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);
	
	std::vector<XYZ> derivatives = ComputeRationalCurveDerivatives(curve, 2, paramT);
	XYZ d1 = derivatives[1];
	XYZ d2 = derivatives[2];
	if (MathUtils::IsAlmostEqualTo(d1.Length(), 1.0))
	{
		return d2.Length();
	}
	double numerator = d1.CrossProduct(d2).Length();
	double denominator = pow(d1.Length(), 3);
	return numerator / denominator;
}

LNLib::XYZ LNLib::NurbsCurve::Normal(const LN_NurbsCurve& curve, CurveNormal normalType, double paramT)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	std::vector<XYZ> derivatives = ComputeRationalCurveDerivatives(curve, 2, paramT);
	XYZ tangent = derivatives[1];
	XYZ der2 = derivatives[2];
	if (MathUtils::IsAlmostEqualTo(tangent.Length(), 1.0))
	{
		XYZ curveNormal = der2 / der2.Length();
		if (normalType == CurveNormal::Normal)
		{
			return curveNormal;
		}
		else
		{
			return tangent.CrossProduct(curveNormal);
		}
	}
	else
	{
		XYZ b = tangent.CrossProduct(der2) / (tangent.CrossProduct(der2).Length());
		if (normalType == CurveNormal::Binormal)
		{
			return b;
		}
		else
		{
			return b.Normalize().CrossProduct(tangent.Normalize());
		}
	}
}

std::vector<LNLib::XYZ> LNLib::NurbsCurve::ProjectNormal(const LN_NurbsCurve& curve)
{
	const std::vector<double>& knotVector = curve.KnotVector;
	int size = knotVector.size();
	int m = size - 1;

	double v0 = knotVector[0];
	double vm = knotVector[m];

	std::vector<XYZ> Blist(size);
	XYZ T0 = ComputeRationalCurveDerivatives(curve, 1, v0)[1].Normalize();

	bool flag = true;
	while (flag)
	{
		bool needReCal = false;
		Blist[0] = XYZ::CreateRandomOrthogonal(T0);

		for (int i = 1; i <= m; i++)
		{
			XYZ Ti = ComputeRationalCurveDerivatives(curve, 1, knotVector[i])[1].Normalize();

			needReCal = Ti.CrossProduct(Blist[i - 1]).IsZero();
			if (needReCal)
			{
				break;
			}

			XYZ bi = Blist[i - 1] - (Blist[i - 1].DotProduct(Ti)) * Ti;
			Blist[i] = bi.Normalize();
		}
		if (!needReCal)
		{
			flag = false;
		}
	}

	bool isClosed = IsClosed(curve);
	if (isClosed)
	{
		Blist[m] = Blist[0];

		std::vector<XYZ> Baver(m);
		Baver[m] = Blist[0];

		for (int i = m - 1; i >= 1; i--)
		{
			XYZ Ti1 = ComputeRationalCurveDerivatives(curve, 1, knotVector[i + 1])[1].Normalize();
			XYZ bi = Baver[i + 1] + Ti1;
			Baver[i] = bi.Normalize();
		}

		for (int i = 1; i <= m - 1; i++)
		{
			Blist[i] = 0.5 * (Blist[i] + Baver[i]);
		}
	}
	return Blist;
}

double LNLib::NurbsCurve::Torsion(const LN_NurbsCurve& curve, double paramT)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	std::vector<XYZ> derivatives = ComputeRationalCurveDerivatives(curve, 3, paramT);
	XYZ tangent = derivatives[1];
	XYZ der2 = derivatives[2];
	XYZ der3 = derivatives[3];
	if (MathUtils::IsAlmostEqualTo(tangent.Length(), 1.0))
	{
		double numerator = (tangent.CrossProduct(der2)).DotProduct(der3);
		double denominator = der2.DotProduct(der2);
		return numerator / denominator;
	}
	else
	{
		double numerator = (tangent.CrossProduct(der2)).DotProduct(der3);
		double denominator = (tangent.CrossProduct(der2).DotProduct(tangent.CrossProduct(der3)));
		return numerator / denominator;
	}
}


int LNLib::NurbsCurve::InsertKnot(const LN_NurbsCurve& curve, double insertKnot, int times, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(times > 0, "times", "Times must be greater than zero.");

	int knotSpanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnot);
	int originMultiplicity = Polynomials::GetKnotMultiplicity(knotVector, insertKnot);

	if (originMultiplicity + times > degree + 1)
	{
		times = degree - 1 - originMultiplicity;
	}
	if (times <= 0)
	{
		return 0;
	}
	
	std::vector<double> insertedKnotVector(knotVector.size() + times);
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

	std::vector<XYZW> updatedControlPoints(controlPoints.size() + times);
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
		if (degree - j - originMultiplicity > 0)
		{
			updatedControlPoints[knotSpanIndex + times - j - originMultiplicity] = temp[degree - j - originMultiplicity];
		}
	}

	for (int i = L + 1; i < knotSpanIndex - originMultiplicity; i++)
	{
		updatedControlPoints[i] = temp[i - L];
	}

	std::sort(insertedKnotVector.begin(), insertedKnotVector.end());
	result.Degree = degree;
	result.KnotVector = insertedKnotVector;
	result.ControlPoints = updatedControlPoints;

	return times;
}

LNLib::XYZ LNLib::NurbsCurve::GetPointOnCurveByCornerCut(const LN_NurbsCurve& curve, double paramT)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(paramT, knotVector[0], knotVector[knotVector.size() - 1]);

	if (MathUtils::IsAlmostEqualTo(paramT, knotVector[0]))
	{
		return controlPoints[0].ToXYZ(true);
	}
	int n = controlPoints.size() - 1;
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

void LNLib::NurbsCurve::RefineKnotVector(const LN_NurbsCurve& curve, const std::vector<double>& insertKnotElements, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(insertKnotElements.size() > 0, "insertKnotElements", "insertKnotElements size must be greater than zero.");

	int n = controlPoints.size() - 1;
	int m = n + degree + 1;
	int r = insertKnotElements.size() - 1;

	int a = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnotElements[0]);
	int b = Polynomials::GetKnotSpanIndex(degree, knotVector, insertKnotElements[r]) + 1;

	std::vector<double> insertedKnotVector(m + r + 2);
	for (int j = 0; j <= a; j++)
	{
		insertedKnotVector[j] = knotVector[j];
	}
	for (int j = b + degree; j <= m; j++)
	{
		insertedKnotVector[j + r + 1] = knotVector[j];
	}

	std::vector<XYZW> updatedControlPoints(n + r + 2);
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
	result.Degree = degree;
	result.KnotVector = insertedKnotVector;
	result.ControlPoints = updatedControlPoints;
}

std::vector<LNLib::LN_NurbsCurve> LNLib::NurbsCurve::DecomposeToBeziers(const LN_NurbsCurve& curve)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	int knotSize = 2 * (degree + 1);
	std::vector<double> bezierKnots(knotSize);
	for (int i = 0; i < knotSize / 2; i++)
	{
		bezierKnots[i] = 0;
	}	
	for (int i = knotSize / 2; i < knotSize; i++)
	{
		bezierKnots[i] = 1;
	}
	int bezierSize = controlPoints.size() - degree;
	std::vector<LNLib::LN_NurbsCurve> beziers(bezierSize);
	for (int i = 0; i < bezierSize; i++)
	{
		beziers[i].Degree = degree;
		beziers[i].KnotVector = bezierKnots;
		beziers[i].ControlPoints = std::vector<XYZW>(degree + 1);
	}

	int n = controlPoints.size() - 1;
	int m = n + degree + 1;

	int a = degree;
	int b = degree + 1;

	int nb = 0;
	for (int i = 0; i <= degree; i++)
	{
		beziers[nb].ControlPoints[i] = controlPoints[i];
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
					beziers[nb].ControlPoints[k] = alpha * beziers[nb].ControlPoints[k] + (1.0 - alpha) * beziers[nb].ControlPoints[k-1];
				}

				if (b < m)
				{
					beziers[nb + 1].ControlPoints[save] = beziers[nb].ControlPoints[degree];
				}
			}

			nb++;
			if (b < m)
			{
				for (int i = degree - multi; i <= degree; i++)
				{
					beziers[nb].ControlPoints[i] = controlPoints[b - degree + i];
				}

				a = b;
				b += 1;
			}
		}
	}
	return beziers;
}

bool LNLib::NurbsCurve::RemoveKnot(const LN_NurbsCurve& curve, double removeKnot, int times, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(removeKnot, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(times > 0, "times", "Times must be greater than zero.");

	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);
	int n = controlPoints.size() - 1;

	int order = degree + 1;
	int s = Polynomials::GetKnotMultiplicity(knotVector, removeKnot);
	int r = Polynomials::GetKnotSpanIndex(degree, knotVector, removeKnot);

	int first = r - degree;
	int last = r - s;

	std::vector<double> restKnotVector = knotVector;
	int m = n + degree + 1;
	for (int k = r + 1; k <= m; k++)
	{
		restKnotVector[k - times] = restKnotVector[k];
	}
	for (int i = 0; i < times; i++)
	{
		restKnotVector.pop_back();
	}

	std::vector<XYZW> updatedControlPoints = controlPoints;
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
			if (MathUtils::IsLessThanOrEqual(controlPoints[i].Distance(alphai * temp[ii + t + 1] + (1.0 - alphai) * temp[ii - 1]), tol))
			{
				remflag = true;
			}
		}

		if (!remflag)
		{
			break;
		}

		i = first;
		j = last;

		while (j - i > t)
		{
			updatedControlPoints[i] = temp[i - off];
			updatedControlPoints[j] = temp[j - off];
			i = i + 1;
			j = j - 1;
		}

		first = first - 1;
		last = last + 1;
	}

	if (t == 0)
	{
		return false;
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
	result.Degree = degree;
	result.KnotVector = restKnotVector;
	result.ControlPoints = updatedControlPoints;
	return true;
}

void LNLib::NurbsCurve::RemoveExcessiveKnots(const LN_NurbsCurve& curve, LN_NurbsCurve& result)
{
	result = curve;
	auto map = KnotVectorUtils::GetInternalKnotMultiplicityMap(curve.KnotVector);
	for (auto it = map.begin(); it != map.end(); it++)
	{
		double u = it->first;
		int count = it->second;

		LN_NurbsCurve tempResult;
		bool canRemove = NurbsCurve::RemoveKnot(result, u, count, tempResult);
		if (!canRemove)
		{
			continue;
		}
		else
		{
			result = tempResult;
		}
	}
}

void LNLib::NurbsCurve::ElevateDegree(const LN_NurbsCurve& curve, int times, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(times > 0, "times", "Times must be greater than zero.");

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

	int moresize = controlPoints.size() + controlPoints.size() * times;
	std::vector<XYZW> updatedControlPoints(moresize,XYZW(Constants::MaxDistance, Constants::MaxDistance, Constants::MaxDistance,1));
	updatedControlPoints[0] = controlPoints[0];

	std::vector<double> updatedKnotVector(moresize + ph + 1,Constants::MaxDistance);
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

	for (int i = updatedControlPoints.size() - 1; i > 0; i--)
	{
		if (MathUtils::IsAlmostEqualTo(updatedControlPoints[i][0], Constants::MaxDistance) &&
			MathUtils::IsAlmostEqualTo(updatedControlPoints[i][1], Constants::MaxDistance) &&
			MathUtils::IsAlmostEqualTo(updatedControlPoints[i][2], Constants::MaxDistance))
		{
			updatedControlPoints.pop_back();
			continue;
		}
		break;
	}
	for (int i = updatedKnotVector.size() -1 ; i >0; i--)
	{
		if (MathUtils::IsAlmostEqualTo(updatedKnotVector[i], Constants::MaxDistance))
		{
			updatedKnotVector.pop_back();
			continue;
		}
		break;
	}
	result.Degree = ph;
	result.KnotVector = updatedKnotVector;
	result.ControlPoints = updatedControlPoints;
}

bool LNLib::NurbsCurve::ReduceDegree(const LN_NurbsCurve& curve, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	double tol = ValidationUtils::ComputeCurveModifyTolerance(controlPoints);

	int size = controlPoints.size();
	bool isBezier = ValidationUtils::IsValidBezier(degree, size);
	if (!isBezier) return false;

	int r = floor((degree - 1) / 2);
	std::vector<XYZW> updatedControlPoints(degree);
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
		double maxU = (degree - std::sqrt(degree)) / (2 * degree);
		error = error * 0.5 * (1 - alpha[r]) * (MathUtils::Binomial(degree, r) * pow(maxU, r) * pow(1 - maxU, r + 1) * (1 - 2 * maxU));
	}
	if (error > tol) return false;

	auto map = KnotVectorUtils::GetKnotMultiplicityMap(knotVector);
	std::vector<double> updatedKnotVector;
	for (auto it = map.begin(); it != map.end(); ++it)
	{
		double u = it->first;
		int count = it->second - 1;
		for (int i = 0; i < count; i++)
		{
			updatedKnotVector.emplace_back(u);
		}
	}
	result.Degree = degree - 1;
	result.KnotVector = updatedKnotVector;
	result.ControlPoints = updatedControlPoints;
	return true;
}

void LNLib::NurbsCurve::EquallyTessellate(const LN_NurbsCurve& curve, std::vector<XYZ>& tessellatedPoints, std::vector<double>& correspondingKnots)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	std::vector<double> uniqueKv = knotVector;
	uniqueKv.erase(unique(uniqueKv.begin(), uniqueKv.end()), uniqueKv.end());
	int size = uniqueKv.size();
	int num = 100;
	for (int i = 0; i < size - 1; i++)
	{
		double currentU = uniqueKv[i];
		double nextU = uniqueKv[i + 1];

		double step = (nextU - currentU) / (num - 1);
		for (int j = 0; j < num; j++)
		{
			double u = currentU + step * j;
			correspondingKnots.emplace_back(u);
			tessellatedPoints.emplace_back(GetPointOnCurve(curve, u));
		}
	}
}

bool LNLib::NurbsCurve::IsClosed(const LN_NurbsCurve& curve)
{
	std::vector<double> knotVector = curve.KnotVector;
	double first = knotVector[0];
	double end = knotVector[knotVector.size() - 1];

	XYZ startPoint = GetPointOnCurve(curve, first);
	XYZ endPoint = GetPointOnCurve(curve, end);

	double distance = startPoint.Distance(endPoint);
	bool isMatch = MathUtils::IsAlmostEqualTo(distance, 0.0);
	if (isMatch) return true;

	std::vector<XYZW> controlPoints = curve.ControlPoints;
	int n = controlPoints.size() - 1;
	XYZ last = controlPoints[n].ToXYZ(true);

	bool flag = false;
	int index = 0;
	for (int i = 0; i < n; i++)
	{
		XYZ current = controlPoints[i].ToXYZ(true);
		if (last.IsAlmostEqualTo(current))
		{
			index = i;
			flag = true;
			break;
		}
	}

	if (!flag) return false;
	if (index == 0) return true;
	for (int i = index; i >= 0; i--)
	{
		XYZ current = controlPoints[i].ToXYZ(true);
		XYZ another = controlPoints[n - index + i].ToXYZ(true);
		if (!another.IsAlmostEqualTo(current))
		{
			return false;
		}
	}
	return true;
}

double LNLib::NurbsCurve::GetParamOnCurve(const LN_NurbsCurve& curve, const XYZ& givenPoint)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	double minValue = Constants::MaxDistance;

	int maxIterations = 10;
	double paramT = Constants::DoubleEpsilon;
	double minParam = knotVector[0];
	double maxParam = knotVector[knotVector.size() - 1];

	std::vector<XYZ> tessellatedPoints;
	std::vector<double> correspondingKnots;
	EquallyTessellate(curve, tessellatedPoints, correspondingKnots);
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

	bool isClosed = IsClosed(curve);
	double a = minParam;
	double b = maxParam;

	int counters = 0;
	while (counters < maxIterations)
	{
		std::vector<XYZ> derivatives = ComputeRationalCurveDerivatives(curve, 2, paramT);
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

		if (!isClosed)
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
		if (condition4 < Constants::DistanceEpsilon) 
		{
			return paramT;
		}

		paramT = temp;
		counters++;
	}
	return paramT;
}

void LNLib::NurbsCurve::CreateTransformed(const LN_NurbsCurve& curve, const Matrix4d& matrix, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	int size = controlPoints.size();
	std::vector<XYZW> transformedControlPoints(size);

	Matrix4d tempMatrix = matrix;
	for (int i = 0; i < size; i++)
	{
		XYZW temp = controlPoints[i];
		transformedControlPoints[i] = tempMatrix.OfWeightedPoint(temp);
	}
	result = curve;
	result.ControlPoints = transformedControlPoints;
}

void LNLib::NurbsCurve::Reparametrize(const LN_NurbsCurve& curve, double alpha, double beta, double gamma, double delta, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(alpha * delta, gamma * beta), "coefficient", "(alpha * delta - gamma * beta) must be greater than zero");

	std::vector<double> updatedKnotVector(knotVector.size());
	for (int i = 0; i < knotVector.size(); i++)
	{
		updatedKnotVector[i] = (alpha * knotVector[i] + beta) / (gamma * knotVector[i] + delta);
	}

	std::vector<XYZW> updatedControlPoints(controlPoints.size());
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

	result.Degree = degree;
	result.KnotVector = updatedKnotVector;
	result.ControlPoints = updatedControlPoints;
}

void LNLib::NurbsCurve::Reparametrize(const LN_NurbsCurve& curve, double min, double max, LN_NurbsCurve& result)
{
	std::vector<double> knotVector = curve.KnotVector;

	if (MathUtils::IsAlmostEqualTo(min, knotVector[0]) && MathUtils::IsAlmostEqualTo(max, knotVector[knotVector.size() - 1]))
	{
		result = curve;
		return;
	}

	std::vector<double> newKnotVector = KnotVectorUtils::Rescale(knotVector, min, max);
	result = curve;
	result.KnotVector = newKnotVector;
}

void LNLib::NurbsCurve::Reverse(const LN_NurbsCurve& curve, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;;

	int size = knotVector.size();
	std::vector<double> reversedKnotVector(size);
	double min = knotVector[0];
	double max = knotVector[size - 1];

	reversedKnotVector[0] = min;
	for (int i = 1; i < size; i++)
	{
		reversedKnotVector[i] = reversedKnotVector[i - 1] + (knotVector[size - i] - knotVector[size - i - 1]);
	}
	std::reverse(controlPoints.begin(), controlPoints.end());

	result.Degree = curve.Degree;
	result.KnotVector = reversedKnotVector;
	result.ControlPoints = controlPoints;
}

bool LNLib::NurbsCurve::SplitAt(const LN_NurbsCurve& curve, double parameter, LN_NurbsCurve& left, LN_NurbsCurve& right)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	if (MathUtils::IsLessThanOrEqual(parameter, knotVector[degree]) ||
		MathUtils::IsGreaterThanOrEqual(parameter, knotVector[knotVector.size() - degree - 1]))
	{
		return false;
	}

	int multi = Polynomials::GetKnotMultiplicity(knotVector, parameter);
	std::vector<double> insert(degree + 1 - multi, parameter);
	left = curve;
	if (insert.size() > 0)
	{
		LN_NurbsCurve temp;
		NurbsCurve::RefineKnotVector(left, insert, temp);
		left = temp;
	}

	int spanIndex = Polynomials::GetKnotSpanIndex(left.Degree, left.KnotVector, parameter) - degree;
	right.Degree = degree;
	int rControlPoints = left.ControlPoints.size() - spanIndex;
	std::vector<XYZW> rightControlPoints(rControlPoints);
	std::vector<double> rightKnotVector(rControlPoints + degree + 1);
	right.KnotVector = rightKnotVector;
	right.ControlPoints = rightControlPoints;
	for (int i = left.ControlPoints.size() - 1, j = rControlPoints - 1; j >= 0; j--, i--)
	{
		right.ControlPoints[j] = left.ControlPoints[i];
	}

	for (int i = left.KnotVector.size() - 1, j = rControlPoints + degree; j >= 0; j--, i--)
	{
		right.KnotVector[j] = left.KnotVector[i];
	}

	left.Degree = degree;
	left.ControlPoints.resize(spanIndex);
	left.KnotVector.resize(spanIndex + degree + 1);

	return true;
}

bool LNLib::NurbsCurve::Segment(const LN_NurbsCurve& curve, double startParameter, double endParameter, LN_NurbsCurve& segment)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	if (MathUtils::IsLessThanOrEqual(startParameter, knotVector[degree]))
	{
		if (MathUtils::IsGreaterThanOrEqual(endParameter, knotVector[knotVector.size() - degree - 1]))
		{
			segment = curve;
			return true;
		}
		LN_NurbsCurve left;
		LN_NurbsCurve right;
		bool result = SplitAt(curve, endParameter, left, right);
		if (result)
		{
			segment = left;
			return true;
		}
	}
	else
	{
		if (MathUtils::IsGreaterThanOrEqual(endParameter, knotVector[knotVector.size() - degree - 1]))
		{
			LN_NurbsCurve left;
			LN_NurbsCurve right;
			bool result = SplitAt(curve, startParameter, left, right);
			if (result)
			{
				segment = right;
				return true;
			}
		}
		else
		{
			LN_NurbsCurve left;
			LN_NurbsCurve right;
			bool result = SplitAt(curve, startParameter, left, right);
			if (result)
			{
				result = SplitAt(right, endParameter, left, right);
				if (result)
				{
					segment = left;
					return true;
				}
			}
		}
	}
	return false;
}

bool LNLib::NurbsCurve::Merge(const LN_NurbsCurve& left, const LN_NurbsCurve& right, LN_NurbsCurve& result)
{
	int degree_L = left.Degree;
	std::vector<double> knotVector_L = left.KnotVector;
	std::vector<XYZW> controlPoints_L = left.ControlPoints;

	int degree_R = right.Degree;
	std::vector<double> knotVector_R = right.KnotVector;
	std::vector<XYZW> controlPoints_R = right.ControlPoints;

	if (!controlPoints_L[controlPoints_L.size() - 1].IsAlmostEqualTo(controlPoints_R[0]))
	{
		return false;
	}

	int degree = std::max(degree_L, degree_R);
	LN_NurbsCurve  tempL;
	Reparametrize(left, 0.0, 1.0, tempL);
	if (degree > degree_L)
	{
		int times = degree - degree_L;
		ElevateDegree(left, times, tempL);
	}

	LN_NurbsCurve tempR;
	Reparametrize(right, 0.0, 1.0, tempR);
	if (degree > degree_R)
	{
		int times = degree - degree_R;
		ElevateDegree(right, times, tempR);
	}

	int l = Polynomials::GetKnotMultiplicity(tempL.KnotVector, tempL.KnotVector[tempL.KnotVector.size() - 1]);
	int r = Polynomials::GetKnotMultiplicity(tempR.KnotVector, tempR.KnotVector[0]);

	if (l != degree + 1 || r != degree + 1)
	{
		return false;
	}

	int size = tempL.ControlPoints.size() + tempR.ControlPoints.size() - 1;
	std::vector<XYZW> controlPoints(size);
	int ksize = size + degree + 1;
	std::vector<double> knotVector(ksize);

	int i;
	for (i = 0; i < tempL.ControlPoints.size(); i++)
	{
		controlPoints[i] = tempL.ControlPoints[i];
	}
		
	for (; i < size; i++)
	{
		controlPoints[i] = tempR.ControlPoints[i - tempL.ControlPoints.size() + 1];
	}

	double kl = tempL.KnotVector[tempL.KnotVector.size() - 1];
	for (i = 0; i < degree + 1; i++)
	{
		knotVector[i] = tempL.KnotVector[0];
		knotVector[ksize - i - 1] = kl + tempR.KnotVector[tempR.KnotVector.size() - 1];
	}
	
	for (int j = 0; j < degree; j++, i++)
	{
		knotVector[i] = kl;
	}

	auto map = KnotVectorUtils::GetInternalKnotMultiplicityMap(tempR.KnotVector);
	for (auto it = map.begin(); it != map.end(); it++)
	{
		for (int j = 0; j < degree; j++, i++)
		{
			knotVector[i] = kl + it->first;
		}
	}	

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = controlPoints;
	return true;
}

void LNLib::NurbsCurve::Offset(const LN_NurbsCurve& curve, double offset, OffsetType type, LN_NurbsCurve& result)
{
	switch (type)
	{
		case OffsetType::TillerAndHanson:
			LNLib::OffsetByTillerAndHanson(curve, offset, result);
			break;

		case OffsetType::PieglAndTiller:
			LNLib::OffsetByPieglAndTiller(curve, offset, result);
		
		default:
			break;
	}
}


void LNLib::NurbsCurve::CreateLine(const XYZ& start, const XYZ& end, LN_NurbsCurve& result)
{
	VALIDATE_ARGUMENT(!start.IsAlmostEqualTo(end), "end", "start must not be equal to end.");

	int degree = 1;
	std::vector<XYZW> controlPoints(2);
	controlPoints[0] = XYZW(start, 1);
	controlPoints[1] = XYZW(end, 1);
	std::vector<double> knotVector(4);
	knotVector[0] = knotVector[1] = 0;
	knotVector[2] = knotVector[3] = 1;

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = controlPoints;
}

void LNLib::NurbsCurve::CreateCubicHermite(const std::vector<XYZ>& throughPoints, const std::vector<XYZ>& tangents, LN_NurbsCurve& curve)
{
	int n = throughPoints.size();
	VALIDATE_ARGUMENT(n > 3, "throughPoints", "ThroughPoints size must be greater than three.");
	VALIDATE_ARGUMENT(n == tangents.size(), "tangents", "Tangents size must be equal to throughPoints size.");
	
	XYZ startPoint = throughPoints[0];
	XYZ endPoint = throughPoints[throughPoints.size() - 1];

	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);

	bool isCyclePoint = startPoint.IsAlmostEqualTo(endPoint);
	XYZ startTangent = tangents[0];
	XYZ endTangent = tangents[tangents.size()-1];
	bool isCycleTangent = startTangent.IsAlmostEqualTo(endTangent);

	int kn = 2 * n;
	int knotSize = kn + 3 + 1;
	std::vector<double> knotVector(knotSize);
	for (int i = 2,  j = 0; i < kn + 2; i += 2, j++)
	{
		knotVector[i] = knotVector[i + 1] = uk[j];
	}
	
	if (isCyclePoint && isCycleTangent)
	{
		knotVector[0] = knotVector[1] = uk[0] - (uk[n - 1] - uk[n - 2]);
		knotVector[kn + 2] = knotVector[kn + 3] = uk[n - 1] + uk[1] - uk[0];
	}
	else if (isCyclePoint && !isCycleTangent)
	{
		knotVector[0] = uk[0] - (uk[n - 1] - uk[n - 2]);
		knotVector[1] = knotVector[2];
		knotVector[kn + 2] = knotVector[kn];
		knotVector[kn + 3] = uk[n - 1] + uk[1] - uk[0];
	}
	else
	{
		knotVector[0] = knotVector[1] = knotVector[2];
		knotVector[kn + 2] = knotVector[kn + 3] = knotVector[kn];
	}

	std::vector<XYZW> controlPoints(kn);
	for (int j = 0, coef = 0; j < kn; j += 2, coef ++)
	{
		double i1 = knotVector[j + 3] - knotVector[j + 1];
		double i2 = knotVector[j + 4] - knotVector[j + 2];

		controlPoints[j] = XYZW((throughPoints[coef] - i1 * tangents[coef] / 3.0),1);
		controlPoints[j + 1] = XYZW((throughPoints[coef] + i2 * tangents[coef] / 3.0), 1);

	}

	curve.Degree = 3;
	curve.KnotVector = knotVector;
	curve.ControlPoints = controlPoints;
}

bool LNLib::NurbsCurve::CreateArc(const XYZ& center, const XYZ& xAxis, const XYZ& yAxis, double startRad, double endRad, double xRadius, double yRadius, LN_NurbsCurve& curve)
{
	VALIDATE_ARGUMENT(!xAxis.IsZero(), "xAxis", "xAxis must not be zero vector.");
	VALIDATE_ARGUMENT(!yAxis.IsZero(), "yAxis", "yAxis must not be zero vector.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(endRad, startRad), "endRad", "endRad must be greater than startRad.");
	double theta = endRad - startRad;
	VALIDATE_ARGUMENT_RANGE(theta, 0, 2 * Constants::Pi);
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(xRadius, 0.0), "xRadius", "xRadius must be greater than zero.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(yRadius, 0.0), "yRadius", "yRadius must be greater than zero.");

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

	int degree = 2;
	std::vector<double> knotVector(n + degree + 1 + 1);
	std::vector<XYZW> controlPoints(n + 1);

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
		knotVector[3] = knotVector[4] = 1.0 / 3;
		knotVector[5] = knotVector[6] = 2.0 / 3;
		break;
	case 4:
		knotVector[3] = knotVector[4] = 0.25;
		knotVector[5] = knotVector[6] = 0.5;
		knotVector[7] = knotVector[8] = 0.75;
		break;
	}
	curve.Degree = degree;
	curve.KnotVector = knotVector;
	curve.ControlPoints = controlPoints;
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
			double a = std::sqrt(alf2 / (1.0 - alf2));
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
			double a = std::sqrt(alf2 / (1.0 - alf2));
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
	insertWeight = std::sqrt(1 + projectPointWeight) / 2.0;
}

bool LNLib::NurbsCurve::CreateOpenConic(const XYZ& start, const XYZ& startTangent, const XYZ& end, const XYZ& endTangent, const XYZ& pointOnConic, LN_NurbsCurve& curve)
{
	VALIDATE_ARGUMENT(!startTangent.IsZero(), "startTangent", "StartTangent must not be zero vector.");
	VALIDATE_ARGUMENT(!endTangent.IsZero(), "endTangent", "EndTangent must not be zero vector.");

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

	int degree = 2;
	std::vector<double> knotVector(j + degree + 1);
	std::vector<XYZW> controlPoints(n + 1);

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
		curve.Degree = degree;
		curve.KnotVector = knotVector;
		curve.ControlPoints = controlPoints;
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
		curve.Degree = degree;
		curve.KnotVector = knotVector;
		curve.ControlPoints = controlPoints;
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
		curve.Degree = degree;
		curve.KnotVector = knotVector;
		curve.ControlPoints = controlPoints;
		return true;
	}
	return false;
}

void LNLib::NurbsCurve::GlobalInterpolation(int degree, const std::vector<XYZ>& throughPoints, LN_NurbsCurve& curve, const std::vector<double>& params)
{
	VALIDATE_ARGUMENT(degree >= 0 && degree <= Constants::NURBSMaxDegree, "degree", "Degree must be greater than or equal zero and not exceed the maximun degree.");
	VALIDATE_ARGUMENT(throughPoints.size() > degree, "throughPoints", "ThroughPoints size must be greater than degree.");
	
	int size = throughPoints.size();
	int n = size - 1;

	std::vector<double> uk(size);
	if (params.size() == 0)
	{
		uk = Interpolation::GetChordParameterization(throughPoints);
	}
	else
	{
		VALIDATE_ARGUMENT(params.size() == size, "params", "Params size must be equal to throughPoints size.");
		uk = params;
	}
	std::vector<double> knotVector = Interpolation::AverageKnotVector(degree, uk);

	std::vector<std::vector<double>> A(size, std::vector<double>(size));
	for (int i = 1; i < n; i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, uk[i]);
		double basis[Constants::NURBSMaxDegree + 1];
		Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i], basis);

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

	std::vector<XYZW> controlPoints(size);
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
	curve.Degree = degree;
	curve.KnotVector = knotVector;
	curve.ControlPoints = controlPoints;
}

void LNLib::NurbsCurve::GlobalInterpolation(int degree, const std::vector<XYZ>& throughPoints, const std::vector<XYZ>& tangents, double tangentFactor, LN_NurbsCurve& curve)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must be greater than zero.");
	VALIDATE_ARGUMENT(throughPoints.size() > degree, "throughPoints", "ThroughPoints size must be greater than degree.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(tangentFactor, 0.0), "tangentFactor", "TangentFactor must be greater than zero.");

	std::vector<XYZ> unitTangents(tangents.size());
	for (int i = 0; i < tangents.size(); i++)
	{
		unitTangents[i] = const_cast<XYZ&>(tangents[i]).Normalize();
	}

	int size = throughPoints.size();
	int n = 2 * size;

	std::vector<XYZW> controlPoints(n);
	std::vector<double> knotVector(n + degree + 1);

	double d = Interpolation::GetTotalChordLength(throughPoints);
	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
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
		double basis[Constants::NURBSMaxDegree + 1];
		Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i], basis);
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
	curve.Degree = degree;
	curve.KnotVector = knotVector;
	curve.ControlPoints = controlPoints;
}

bool LNLib::NurbsCurve::CubicLocalInterpolation(const std::vector<XYZ>& throughPoints, LN_NurbsCurve& curve)
{
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints size must be greater than zero.");

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

		double alpha = (-b + std::sqrt(b * b - 4 * a * c)) / (2 * a);

		XYZ pk0 = throughPoints[k];
		XYZ pk1 = p0 + (1 / 3) * alpha * t0;
		XYZ pk2 = p3 - (1 / 3) * alpha * t3;

		uk[k + 1] = uk[k] + 3 * (pk1 - pk0).Length();

		tempControlPoints.emplace_back(XYZW(pk1, 1));
		tempControlPoints.emplace_back(XYZW(pk2, 1));
	}

	int kvSize = 2 * (degree + 1) + 2 * (size - 1);
	std::vector<double> knotVector(kvSize);
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
	std::vector<XYZW> controlPoints(tSize + 2);
	controlPoints[0] = XYZW(throughPoints[0], 1);
	controlPoints[tSize + 2 - 1] = XYZW(throughPoints[n], 1);
	for (int i = 0; i < tSize; i++)
	{
		controlPoints[i + 1] = tempControlPoints[i];
	}
	curve.Degree = degree;
	curve.KnotVector = knotVector;
	curve.ControlPoints = controlPoints;
	return true;
}

bool LNLib::NurbsCurve::LeastSquaresApproximation(int degree, const std::vector<XYZ>& throughPoints, int controlPointsCount, LN_NurbsCurve& curve)
{
	VALIDATE_ARGUMENT(degree >= 0 && degree <= Constants::NURBSMaxDegree, "degree", "Degree must be greater than or equal zero and not exceed the maximun degree.");
	VALIDATE_ARGUMENT(controlPointsCount > 0, "controlPointsCount", "controlPointsCount must be greater than zero.");

	int n = controlPointsCount;
	int m = throughPoints.size();
	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	std::vector<double> knotVector(n + degree + 1, 1.0);
	
	double d = (double)m / (double)n;
	for (int j = 0; j <= degree; j++)
	{
		knotVector[j] = 0.0;
	}
	for (int j = 1; j < n - degree; j++)
	{
		knotVector[degree + j] = 0.0;
		for (int k = j; k < j + degree; k++)
		{
			int i1 = k * d;
			double a = k * d - i1;
			int i2 = ((k - 1) * d);
			knotVector[degree + j] += a * uk[i2] + (1 - a) * uk[i1];
		}
		knotVector[degree + j] /= degree;
	}

	std::vector<XYZW> controlPoints(n);
	std::vector<XYZ> R(n);
	std::vector<XYZ> rk(m);
	std::vector<double> funs(degree + 1);
	std::vector<std::vector<double>> N(m, std::vector<double>(n));
	R[0] = throughPoints[0];
	R[n - 1] = throughPoints[m - 1];
	N[0][0] = 1.0;
	N[m - 1][n - 1] = 1.0;

	for (int i = 0; i < m; i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, uk[i]);
		double basis[Constants::NURBSMaxDegree + 1];
		Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i], basis);
		for (int j = 0; j <= degree; j++)
		{
			N[i][spanIndex - degree + j] = basis[j];
		}
		rk[i] = throughPoints[i] - N[i][0] * throughPoints[0] - N[i][n - 1] * throughPoints[m - 1];
	}

	for (int i = 0; i < n; i++) 
	{
		R[i] = XYZ();
		for (int j = 0; j < m; j++)
		{
			R[i] += N[j][i] * rk[j];
		}
			
		if (R[i].IsAlmostEqualTo(XYZ()))
		{
			return false;
		}
	}
	
	if (n - 2 > 0) 
	{
		std::vector<std::vector<double>> X(n - 2, std::vector<double>(3));
		std::vector<std::vector<double>> B(n - 2, std::vector<double>(3));
		std::vector<std::vector<double>> Ns(m - 2, std::vector<double>(n-2));
		for (int i = 0; i < n - 2; i++) 
		{
			for (int j = 0; j < 3; j++)
			{
				B[i][j] = R[i + 1][j];
			}	
		}
		for (int i = 1; i <= m - 2; i++)
		{
			for (int j = 1; j <= n - 2; j++)
			{
				Ns[i - 1][j - 1] = N[i][j];
			}
		}
		std::vector<std::vector<double>> NsT;
		MathUtils::Transpose(Ns, NsT);
		auto NsTNs = MathUtils::MatrixMultiply(NsT, Ns);
		X = MathUtils::SolveLinearSystem(NsTNs, B);

		for (int i = 0; i < n - 2; i++)
		{
			double x = X[i][0];
			double y = X[i][1];
			double z = X[i][2];
			controlPoints[i+1] = XYZW(XYZ(x,y,z),1);
		}
	}
	controlPoints[0] = XYZW(throughPoints[0],1);
	controlPoints[n - 1] = XYZW(throughPoints[m - 1],1);

	curve.Degree = degree;
	curve.KnotVector = knotVector;
	curve.ControlPoints = controlPoints;
	return true;
}

bool LNLib::NurbsCurve::WeightedAndContrainedLeastSquaresApproximation(int degree, const std::vector<XYZ>& throughPoints, const std::vector<double>& throughPointWeights, const std::vector<XYZ>& tangents, const std::vector<int>& tangentIndices, const std::vector<double>& tangentWeights, int controlPointsCount, LN_NurbsCurve& curve)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must be greater than zero.");
	int size = throughPoints.size();
	VALIDATE_ARGUMENT(size > degree, "throughPoints", "ThroughPoints size must be greater than degree.");
	VALIDATE_ARGUMENT(throughPointWeights.size() == size, "weights", "Weights size must be equal to throughPoints size.");

	
	int ru = -1;
	int rc = -1;
	int r = size - 1;
	for (int i = 0; i <= r; i++)
	{
		if (MathUtils::IsGreaterThan(throughPointWeights[i], 0.0))
		{
			ru++;
		}
		else
		{
			rc++;
		}
	}
	int su = -1;
	int sc = -1;
	int s = tangents.size() - 1;
	for (int j = 0; j <= s; j++)
	{
		if (MathUtils::IsGreaterThan(tangentWeights[j], 0.0))
		{
			su++;
		}
		else
		{
			sc++;
		}
	}
	int mu = ru + su + 1;
	int mc = rc + sc + 1;

	int n = controlPointsCount - 1;
	if (mc >= n || mc + n >= mu + 1)
	{
		return false;
	}

	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	std::vector<double> knotVector = Interpolation::ComputeKnotVector(degree, size, controlPointsCount, uk);
	std::vector<XYZW> controlPoints(controlPointsCount);

	int j = 0;
	int mu2 = 0;
	int mc2 = 0;

	std::vector<std::vector<double>> N(mu + 1, std::vector<double>(n + 1));
	std::vector<std::vector<double>> M(mc + 1, std::vector<double>(n + 1));
	std::vector<XYZ> S(mu + 1);
	std::vector<XYZ> T(mc + 1);

	std::vector<double> W(mu + 1);
	std::vector<std::vector<double>> funs(2, std::vector<double>(degree + 1));

	for (int i = 0; i <= r; i++)
	{
		int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, uk[i]);
		
		bool dflag = false;
		if (j <= s)
		{
			if (i == tangentIndices[j])
			{
				dflag = true;
			}
		}
		if (!dflag)
		{
			std::vector<double> temp(degree + 1);
			Polynomials::BasisFunctions(spanIndex, degree, knotVector, uk[i], temp.data());
			funs[0] = temp;
		}
		else
		{
			funs = Polynomials::BasisFunctionsDerivatives(spanIndex, degree, 1, knotVector, uk[i]);
		}
		if (MathUtils::IsGreaterThan(throughPointWeights[i], 0.0))
		{
			W[mu2] = throughPointWeights[i];
			//
			// My question: 
			// 
			// N is (mu + 1) rows, (n + 1) columns;
			// funs is 2 rows, (degree + 1) columns;
			// How to match column data? Is controlPointsCount must equal degree + 1?
			//
			for (int z = 0; z < funs[0].size(); z++)
			{
				N[mu2][z] = funs[0][z];
			}
			S[mu2] = W[mu2] * throughPoints[i];
			mu2++;
		}
		else
		{
			for (int z = 0; z < funs[0].size(); z++)
			{
				M[mc2][z] = funs[0][z];
			}
			T[mc2] = throughPoints[i];
			mc2++;
		}
		if (dflag)
		{
			if (MathUtils::IsGreaterThan(tangentWeights[j], 0.0))
			{
				W[mu2] = tangentWeights[j];
				for (int z = 0; z < funs[1].size(); z++)
				{
					N[mu2][z] = funs[1][z];
				}
				S[mu2] = W[mu2] * tangents[j];
				mu2++;
			}
			else
			{
				for (int z = 0; z < funs[1].size(); z++)
				{
					M[mc2][z] = funs[1][z];
				}
				T[mc2] = tangents[j];
				mc2++;
			}
			j++;
		}
	}
	
	std::vector<std::vector<double>> tN;
	MathUtils::Transpose(N, tN);

	std::vector<std::vector<double>> W_ = MathUtils::MakeDiagonal(W.size());
	std::vector<std::vector<double>> tNW = MathUtils::MatrixMultiply(tN, W_);
	std::vector<std::vector<double>> tNWN = MathUtils::MatrixMultiply(tNW, N);

	std::vector<std::vector<double>> S_(S.size(), std::vector<double>(3));
	for (int i = 0; i < S.size(); i++)
	{
		for (int j = 0; j < 3; j++)
		{
			S_[i][j] = S[i][j];
		}
	}
	std::vector<std::vector<double>> tNWS = MathUtils::MatrixMultiply(tNW, S_);
	
	if (mc < 0)
	{
		std::vector<std::vector<double>> result = MathUtils::SolveLinearSystem(tNWN, tNWS);
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
	else
	{
		std::vector<std::vector<double>> inv_tNWN;
		bool canInverse = MathUtils::MakeInverse(tNWN, inv_tNWN);
		if (!canInverse) return false;

		std::vector<std::vector<double>> tM;
		MathUtils::Transpose(M, tM);

		std::vector<std::vector<double>> Minv_tNWN = MathUtils::MatrixMultiply(M, inv_tNWN);
		std::vector<std::vector<double>> Minv_tNWN_tM = MathUtils::MatrixMultiply(Minv_tNWN, tM);

		std::vector<std::vector<double>> Minv_tNWN_tNWS = MathUtils::MatrixMultiply(Minv_tNWN, tNWS);
		std::vector<std::vector<double>> Minv_tNWN_tNWS_T(T.size(), std::vector<double>(3));
		for (int i = 0; i < T.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				Minv_tNWN_tNWS_T[i][j] = Minv_tNWN_tNWS[i][j] - T[i][j];
			}
		}

		std::vector<std::vector<double>> A = MathUtils::SolveLinearSystem(Minv_tNWN_tM, Minv_tNWN_tNWS_T);
		std::vector<std::vector<double>> tMA = MathUtils::MatrixMultiply(tM, A);
		std::vector<std::vector<double>> tNWS_tMA(tNWS.size(), std::vector<double>(3));
		for (int i = 0; i < tNWS.size(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				tNWS_tMA[i][j] = tNWS[i][j] - tMA[i][j];
			}
		}
		std::vector<std::vector<double>> result = MathUtils::MatrixMultiply(inv_tNWN, tNWS_tMA);
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

	curve.Degree = degree;
	curve.KnotVector = knotVector;
	curve.ControlPoints = controlPoints;

	return true;
}

double LNLib::NurbsCurve::ComputerRemoveKnotErrorBound(const LN_NurbsCurve& curve, int removalIndex)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(removalIndex, 0, knotVector.size()-1);

	int ord = degree + 1;
	int r = removalIndex;
	double u = knotVector[r];
	int s = Polynomials::GetKnotMultiplicity(knotVector, u);
	int last = r - s;
	int first = r - degree;
	int off = first - 1;

	std::vector<XYZW> temp(knotVector.size());
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

void LNLib::NurbsCurve::RemoveKnotsByGivenBound(const LN_NurbsCurve& curve, const std::vector<double> params, std::vector<double>& errors, double maxError, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(params.size() > 0, "params", "Params size must be greater than zero.");
	VALIDATE_ARGUMENT(params.size() == errors.size(), "errors", "Errors size must be equal to params size.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(maxError,0.0), "maxError", "Maxerror must be greater than zero.");

	int knotSize = knotVector.size();
	std::vector<double> Br(knotSize, Constants::MaxDistance);
	std::vector<int> S(knotSize, 0);
	std::vector<int> Nl(knotSize);
	std::vector<int> Nr(knotSize, params.size() - 1);

	std::vector<double> uk = params;
	int ukSize = uk.size();
	std::vector<double> NewError(ukSize);
	std::vector<double> temp(ukSize);

	int i, BrMinIndex = 0;
	int r, s, Rstart, Rend, k;
	double BrMin, u;

	int controlPointsSize = controlPoints.size();
	int n = controlPointsSize - 1;
	s = 1;
	for (int i = degree + 1; i < controlPointsSize; i++)
	{
		if (MathUtils::IsLessThan(knotVector[i],knotVector[i + 1]))
		{
			Br[i] = ComputerRemoveKnotErrorBound(curve, i);
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

	std::vector<double> updatedKnotVector;
	std::vector<XYZW> updatedControlPoints;

	while (1)
	{
		double minStandard = Constants::MaxDistance;
		for (int i = 0; i < Br.size(); i++)
		{
			if (MathUtils::IsLessThan(Br[i],minStandard))
			{
				BrMinIndex = i;
				minStandard = Br[i];
			}
		}
		BrMin = Br[BrMinIndex];

		if (BrMin == Constants::MaxDistance)
		{
			updatedKnotVector = tempU;
			updatedControlPoints = tempCP;
			break;
		}
		r = BrMinIndex;
		s = S[BrMinIndex];

		Rstart = std::max(r - degree, degree + 1);
		Rend = std::min(r + degree - S[r + degree] + 1, n);
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
			temp[i] = NewError[i] + errors[i];
			if (MathUtils::IsGreaterThan(temp[i],maxError))
			{
				removable = false;
				Br[r] = Constants::MaxDistance;
				updatedKnotVector = tempU;
				updatedControlPoints = tempCP;
				break;
			}
		}

		if (removable)
		{
			LN_NurbsCurve tc;
			tc.Degree = degree;
			tc.KnotVector = tempU;
			tc.ControlPoints = tempCP;

			LN_NurbsCurve newtc;
			RemoveKnot(tc, tempU[r], 1, newtc);
			std::vector<double> tempNewU = newtc.KnotVector;
			std::vector<XYZW> tempNewCP = newtc.ControlPoints;

			controlPointsSize = tempNewCP.size();
			for (int i = Rstart; i <= Rend; i++)
			{
				errors[i] = temp[i];
			}

			if (controlPointsSize <= degree + 1)
			{
				updatedKnotVector = tempU;
				updatedControlPoints = tempCP;
				break;
			}

			Rstart = Nl[r - degree - 1];
			Rend = Nr[r - S[r]];
			int spanIndex = 0;
			int oldspanIndex = -1;
			for (int k = Rstart; k <= Rend; k++)
			{
				spanIndex = Polynomials::GetKnotSpanIndex(degree, tempNewU, uk[k]);
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
			Rend = std::min(r + degree - S[r] + 1, controlPointsSize);
			s = S[Rstart];
			for (i = Rstart; i <= Rend; i++)
			{
				if (MathUtils::IsLessThan(tempNewU[i],tempNewU[i + 1]))
				{
					Br[i] = ComputerRemoveKnotErrorBound(newtc, i);
					S[i] = s;
					s = 1;
				}
				else {
					Br[i] = Constants::MaxDistance;
					S[i] = 1;
					s++;
				}
			}
			for (int i = Rend + 1; i < Br.size() - 1; i++)
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

	result.Degree = degree;
	result.KnotVector = updatedKnotVector;
	result.ControlPoints = updatedControlPoints;
}

void LNLib::NurbsCurve::GlobalApproximationByErrorBound(int degree, const std::vector<XYZ>& throughPoints, double maxError, LN_NurbsCurve& result)
{
	VALIDATE_ARGUMENT(degree > 0, "degree", "Degree must be greater than zero.");
	VALIDATE_ARGUMENT(throughPoints.size() > degree, "throughPoints", "ThroughPoints size must be greater than degree.");
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(maxError, 0.0), "maxError", "Maxerror must be greater than zero.");

	std::vector<double> uk = Interpolation::GetChordParameterization(throughPoints);
	int size = throughPoints.size();
	std::vector<double> errors(size, 0.0);

	std::vector<XYZW> controlPoints(size);
	int m = (size-1) + 1 + 1;
	std::vector<double> knotVector(m + 1);

	for (int i = 0; i < size; i++) 
	{
		knotVector[i + 1] = uk[i];
	}
	knotVector[0] = 0;
	knotVector[knotVector.size() - 1] = 1.0;
	
	for (int i = 0; i < size; i++)
	{
		controlPoints[i] = XYZW(throughPoints[i],1);
	}

	LN_NurbsCurve tc;
	tc.Degree = 1;
	tc.KnotVector = knotVector;
	tc.ControlPoints = controlPoints;

	LN_NurbsCurve newtc;
	ElevateDegree(tc, degree - 1, newtc);
	RemoveKnotsByGivenBound(newtc, uk, errors, maxError, result);
}

bool LNLib::NurbsCurve::FitWithConic(const std::vector<XYZ>& throughPoints, int startPointIndex, int endPointIndex, const XYZ& startTangent, const XYZ& endTangent, double maxError, std::vector<XYZW>& middleControlPoints)
{
	VALIDATE_ARGUMENT(throughPoints.size() > 0, "throughPoints", "ThroughPoints size must be greater than zero.");
	VALIDATE_ARGUMENT_RANGE(startPointIndex, 0, throughPoints.size() - 1);
	VALIDATE_ARGUMENT_RANGE(endPointIndex, startPointIndex + 1, throughPoints.size() - 1);
	VALIDATE_ARGUMENT(!startTangent.IsZero(), "startTangent", "StartTangent must not be zero vector.");
	VALIDATE_ARGUMENT(!endTangent.IsZero(), "endTangent", "EndTangent must not be zero vector.");

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
	else if (type == CurveCurveIntersectionType::Skew || type == CurveCurveIntersectionType::Parallel)
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
	XYZ dummy;
	for (int i = startPointIndex + 1; i <= endPointIndex - 1; i++)
	{
		XYZ V1 = throughPoints[i] - R;
		type = Intersection::ComputeRays(startPoint, V, R, V1, alf1, alf2, dummy);
		if (type != CurveCurveIntersectionType::Intersecting || 
			MathUtils::IsLessThanOrEqual(alf1, 0.0) ||
			MathUtils::IsGreaterThanOrEqual(alf1, 1.0) ||
			MathUtils::IsLessThanOrEqual(alf2, 0.0))
		{
			return false;	
		}
		double wi = 0.0;
		if (CreateOneConicArc(startPoint, V, R, V1, throughPoints[i], dummy, wi))
		{
			s = s + wi / (1 + wi);
		}
	}
	
	s = s / (endPointIndex - startPointIndex - 1);
	double w = s / (1.0 - s);

	std::vector<XYZW> controlPoints = { XYZW(startPoint,1), XYZW(R,w), XYZW(endPoint,1) };
	std::vector<double> knotVectors = { 0,0,0,1.0,1.0,1.0 };
	int degree = 2;
	for (int k = startPointIndex + 1; k <= endPointIndex - 1; k++)
	{
		XYZ tp = throughPoints[k];
		LN_NurbsCurve tc;
		tc.Degree = degree;
		tc.KnotVector = knotVectors;
		tc.ControlPoints = controlPoints;

		double param = GetParamOnCurve(tc, tp);
		XYZ point = GetPointOnCurve(tc, param);
		if (MathUtils::IsGreaterThan(tp.Distance(point), maxError))
		{
			return false;
		}
	}
	middleControlPoints.emplace_back(XYZW(R, w));
	return true;
}

bool LNLib::NurbsCurve::FitWithCubic(const std::vector<XYZ>& throughPoints, int startPointIndex, int endPointIndex, const XYZ& startTangent, const XYZ& endTangent, double maxError, std::vector<XYZW>& middleControlPoints)
{
	VALIDATE_ARGUMENT(throughPoints.size() >= 3, "throughPoints", "ThroughPoints size must be greater than 2.");
	VALIDATE_ARGUMENT_RANGE(startPointIndex, 0, throughPoints.size() - 1);
	VALIDATE_ARGUMENT_RANGE(endPointIndex, startPointIndex + 1, throughPoints.size() - 1);
	VALIDATE_ARGUMENT(!startTangent.IsZero(), "startTangent", "StartTangent must not be zero vector.");
	VALIDATE_ARGUMENT(!endTangent.IsZero(), "endTangent", "EndTangent must not be zero vector.");

	int size = throughPoints.size();
	XYZ startPoint = throughPoints[startPointIndex];
	XYZ endPoint = throughPoints[endPointIndex];
	double alpha = 0.0;
	double beta = 0.0;
	if (endPointIndex - startPointIndex == 1)
	{
		XYZ Dks(0, 0, 0);
		XYZ Dke(0, 0, 0);
		std::vector<XYZ> tangents = Interpolation::ComputeTangent(throughPoints);
		if (startPointIndex == 0)
		{
			Dks = tangents[startPointIndex];
		}
		else
		{
			Dks = (endPoint.Distance(startPoint) / startPoint.Distance(throughPoints[startPointIndex - 1])) * startTangent;
		}
		if (endPointIndex == size - 1)
		{
			Dke = tangents[endPointIndex];
		}
		else
		{
			Dke = (throughPoints[endPointIndex + 1].Distance(endPoint) / endPoint.Distance(startPoint)) * endTangent;
		}
		alpha = Dks.Length() / 3.0;
		beta = -Dke.Length() / 3.0;
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
		XYZW P1((2 * startPoint + endPoint) / 3.0, 1);
		XYZW P2((startPoint + 2 * endPoint) / 3.0, 1);
		middleControlPoints.emplace_back(P1);
		middleControlPoints.emplace_back(P2);
		return true;
	}
	std::vector<XYZ> newThroughPoints(endPointIndex - startPointIndex + 1);
	for (int i = startPointIndex; i <= endPointIndex; i++)
	{
		newThroughPoints.emplace_back(throughPoints[i]);
	}
	std::vector<double> uh = Interpolation::GetChordParameterization(newThroughPoints);
	std::vector<double> alphak(dk);
	std::vector<double> betak(dk);
	for (int k = 1; k < dk; k++)
	{
		XYZ normalPi = (endPoint - startPoint).Normalize().CrossProduct(startTangent);
		XYZ intersectPoint;
		auto type = Intersection::ComputeLineAndPlane(normalPi, startPoint, throughPoints[k + startPointIndex], endTangent, intersectPoint);
		if (type == LinePlaneIntersectionType::On)
		{
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
				alphak[k] = ak;
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
			if (!(type0 == LinePlaneIntersectionType::Intersecting)) continue;
			double param0 = 0.0; double param1 = 0.0;
			XYZ Pc(0, 0, 0);
			CurveCurveIntersectionType type1 = Intersection::ComputeRays(startPoint, endPoint - startPoint, Pd, startTangent, param0, param1, Pc);
			if (!(type1 == CurveCurveIntersectionType::Intersecting)) continue;
			double gamma = Pc.Distance(endPoint) / startPoint.Distance(endPoint);
			if (MathUtils::IsLessThan(gamma, 0.0) || MathUtils::IsGreaterThan(gamma, 1.0))
			{
				return false;
			}
			uh[k] = MathUtils::ComputerCubicEquationsWithOneVariable(2, -3, 0, 1 - gamma);
			if (MathUtils::IsLessThan(uh[k], 0.0) || MathUtils::IsGreaterThan(uh[k], 1.0))
			{
				return false;
			}
			else
			{
				double a = Pc.Distance(Pd);
				double b = -Pd.Distance(throughPoints[k + startPointIndex]);
				double b13 = Polynomials::Bernstein(1, 3, uh[k]);
				double b23 = Polynomials::Bernstein(2, 3, uh[k]);
				alphak[k] = a / b13;
				betak[k] = b / b23;
			}
		}
	}
	alpha = 0.0;
	beta = 0.0;
	for (int k = 1; k < dk; k++)
	{
		alpha += alphak[k];
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
		double deltaAk = alphak[k] - alpha;
		double deltaBk = betak[k] - beta;
		double e = (deltaAk * Polynomials::Bernstein(1, 3, u) * startTangent).Distance(deltaBk * Polynomials::Bernstein(2, 3, u) * endTangent);
		if (MathUtils::IsLessThanOrEqual(e, maxError)) continue;
		else
		{
			std::vector<XYZ> cps = { XYZ(startPoint), XYZ(P1), XYZ(P2),  XYZ(endPoint) };

			LN_BezierCurve<XYZ> bezierCurve;
			bezierCurve.Degree = 3;
			bezierCurve.ControlPoints = cps;

			XYZ p = BezierCurve::GetPointOnCurveByBernstein(bezierCurve, u);
			double ek = throughPoints[startPointIndex + k].Distance(p);
			if (MathUtils::IsGreaterThan(ek, maxError))
			{
				return false;
			}
		}
	}
	return true;
}

bool LNLib::NurbsCurve::ControlPointReposition(const LN_NurbsCurve& curve, double parameter, int moveIndex, XYZ moveDirection, double moveDistance, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(parameter, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(moveIndex, 0, controlPoints.size() - 1);
	VALIDATE_ARGUMENT(!moveDirection.IsZero(), "moveDirection", "MoveDirection must not be zero vector.");
	VALIDATE_ARGUMENT(!MathUtils::IsAlmostEqualTo(moveDistance,0.0), "moveDistance", "MoveDistance must not be zero.")

	int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, parameter);
	double basis[Constants::NURBSMaxDegree + 1];
	Polynomials::BasisFunctions(spanIndex, degree, knotVector, parameter, basis);
	double Rkp = basis[0];
	if (MathUtils::IsLessThan(Rkp, 0.0))
	{
		return false;
	}
	std::vector<XYZW> updatedControlPoints = controlPoints;
	XYZ movePoint = updatedControlPoints[moveIndex].ToXYZ(true);
	double alpha = moveDistance / (moveDirection.Length() * Rkp);
	XYZ newPoint = movePoint + alpha * moveDirection;
	updatedControlPoints[moveIndex] = XYZW(newPoint, controlPoints[moveIndex].GetW());

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = updatedControlPoints;
	return true;
}

void LNLib::NurbsCurve::WeightModification(const LN_NurbsCurve& curve, double parameter, int moveIndex, double moveDistance, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(parameter, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(moveIndex, 0, controlPoints.size() - 1);
	VALIDATE_ARGUMENT(!MathUtils::IsAlmostEqualTo(moveDistance, 0.0), "moveDistance", "MoveDistance must not be zero.");

	XYZ point = GetPointOnCurve(curve, parameter);
	XYZ movePoint = const_cast<XYZW&>(controlPoints[moveIndex]).ToXYZ(true);
	double distance =  point.Distance(movePoint);
	int spanIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, parameter);
	double basis[Constants::NURBSMaxDegree + 1];
	Polynomials::BasisFunctions(spanIndex, degree, knotVector, parameter, basis);
	double Rkp = basis[0];
	double coefficient = 1 + moveDistance / (Rkp * (distance - moveDistance));
	std::vector<XYZW> updatedControlPoints = controlPoints;
	updatedControlPoints[moveIndex] = XYZW(movePoint, updatedControlPoints[moveIndex].GetW() * coefficient);

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = updatedControlPoints;
}

bool LNLib::NurbsCurve::NeighborWeightsModification(const LN_NurbsCurve& curve, double parameter, int moveIndex, double moveDistance, double scale, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(parameter, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(moveIndex, 0, controlPoints.size() - 1);
	VALIDATE_ARGUMENT(!MathUtils::IsAlmostEqualTo(moveDistance, 0.0), "moveDistance", "MoveDistance must not be zero.");
	VALIDATE_ARGUMENT(!MathUtils::IsAlmostEqualTo(scale, 0.0), "scale", "Scale must not be zero.");

	std::vector<XYZW> tempControlPoints = controlPoints;
	XYZ movePoint1 = const_cast<XYZW&>(tempControlPoints[moveIndex]).ToXYZ(true);
	tempControlPoints[moveIndex] = XYZW(movePoint1, 0.0);
	XYZ movePoint2 = const_cast<XYZW&>(tempControlPoints[moveIndex+1]).ToXYZ(true);
	tempControlPoints[moveIndex+1] = XYZW(movePoint2, 0.0);

	LN_NurbsCurve tc;
	tc.Degree = degree;
	tc.KnotVector = knotVector;
	tc.ControlPoints = tempControlPoints;

	XYZ R = GetPointOnCurve(tc, parameter);
	XYZ controlLeg = movePoint1 - movePoint2;
	double controlLegLength = movePoint2.Distance(movePoint1);

	XYZ P = GetPointOnCurve(curve, parameter);
	XYZ direction = R - P;
	double param0,param1;
	XYZ Q;
	auto type = Intersection::ComputeRays(movePoint1, controlLeg, R, direction, param0, param1, Q);
	if (type != CurveCurveIntersectionType::Intersecting)
	{
		return false;
	}
	XYZ pkq = Q - movePoint1;
	XYZ pk1q = Q - movePoint2;

	double RQ = Q.Distance(R);
	double RP = P.Distance(R);
	direction = (P - R) / RP;
	
	XYZ target = P + moveDistance * direction;
	double Rtarget = RP + moveDistance;
	double qRP = RP / RQ;
	double qRtarget = Rtarget / RQ;

	XYZ A = movePoint1 + qRP * pkq;
	XYZ B = movePoint2 + qRP * pk1q;
	XYZ C = movePoint1 + qRtarget * pkq;
	XYZ D = movePoint2 + qRtarget * pk1q;

	double ak = B.Distance(movePoint2) / controlLegLength;
	double ak1 = A.Distance(movePoint1) / controlLegLength;
	double abk = D.Distance(movePoint2) / controlLegLength;
	double abk1 = C.Distance(movePoint1) / controlLegLength;

	if (MathUtils::IsLessThan(abs(ak), 0.0) || MathUtils::IsLessThan(abs(ak1), 0.0) ||
		MathUtils::IsLessThan(abs(abk), 0.0) || MathUtils::IsLessThan(abs(abk1), 0.0))
	{
		return false;
	}

	double alpha = 1.0 - ak - ak1;
	double beta = 1.0 - abk - abk1;

	double betak = (alpha / ak) / (beta / abk);
	double betak1 = (alpha / ak1) / (beta / abk1);

	std::vector<XYZW> updatedControlPoints = controlPoints;
	updatedControlPoints[moveIndex] = XYZW(movePoint1, updatedControlPoints[moveIndex].GetW() * betak);
	updatedControlPoints[moveIndex + 1] = XYZW(movePoint2, updatedControlPoints[moveIndex + 1].GetW() * betak1);

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = updatedControlPoints;
	return true;
}

void LNLib::NurbsCurve::Warping(const LN_NurbsCurve& curve, const std::vector<double>& warpShape, double warpDistance, const XYZ& planeNormal, double startParameter, double endParameter, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(controlPoints.size() == warpShape.size(), "warpShape", "WarpShape size must be equal to control points size.");
	VALIDATE_ARGUMENT(!MathUtils::IsAlmostEqualTo(warpDistance, 0.0), "warpDistance", "WarpDistance must not be zero.");
	VALIDATE_ARGUMENT(!planeNormal.IsZero(), "planeNormal", "PlaneNormal must not be zero vector.");
	VALIDATE_ARGUMENT_RANGE(startParameter, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(endParameter, startParameter, knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(endParameter, startParameter), "endParameter", "EndParameter must be greater than startParamter.");

	double halfParameter = 0.5 * (startParameter + endParameter);
	XYZ tangent = ComputeRationalCurveDerivatives(curve, 1, halfParameter)[1];
	XYZ normal = tangent.CrossProduct(planeNormal);
	XYZ W = MathUtils::IsGreaterThan(warpDistance,0.0)? normal : -normal;
	std::vector<XYZW> temp = controlPoints;
	std::vector<XYZW> resultControlPoints;
	for (int i = 0; i < temp.size(); i++)
	{
		double weight = temp[i].GetW();
		XYZ currentPoint = temp[i].ToXYZ(true);
		XYZ newPoint = currentPoint + warpShape[i] * abs(warpDistance) * W.Normalize();
		resultControlPoints.emplace_back(XYZW(newPoint, weight));
	}
	
	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = resultControlPoints;
}

bool LNLib::NurbsCurve::Flattening(const LN_NurbsCurve& curve, XYZ lineStartPoint, XYZ lineEndPoint, double startParameter, double endParameter, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(!lineStartPoint.IsAlmostEqualTo(lineEndPoint), "lineEndPoint", "lineEndPoint must not be equal to lineStartPoint.");
	VALIDATE_ARGUMENT_RANGE(startParameter, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(endParameter, startParameter, knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(endParameter, startParameter), "endParameter", "EndParameter must be greater than startParamter.");

	int spanMinIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, startParameter);
	int spanMaxIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, endParameter);
	
	std::map<int, XYZ> selectedControlPoints;
	for (int i = spanMinIndex; i <= spanMaxIndex - degree - 1; i++)
	{
		XYZ p = const_cast<XYZW&>(controlPoints[i]).ToXYZ(true);
		selectedControlPoints.insert({ i, p });
	}
	
	int projectCount = 0;
	std::vector<XYZW> updatedControlPoints = controlPoints;
	XYZ lineVector = (lineEndPoint - lineStartPoint).Normalize();
	double lineLength = lineEndPoint.Distance(lineStartPoint);
	for (auto it = selectedControlPoints.begin(); it != selectedControlPoints.end(); ++it)
	{
		int index = it->first;
		XYZ current = it->second;
		XYZ project;
		bool isProjection = Projection::PointToLine(lineStartPoint, lineEndPoint, current, project);
		if (isProjection)
		{
			projectCount++;
			updatedControlPoints[index] = XYZW(project, updatedControlPoints[index].GetW());
		}
	}

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = updatedControlPoints;

	return projectCount >= degree + 1;
}

void LNLib::NurbsCurve::Bending(const LN_NurbsCurve& curve, double startParameter, double endParameter, XYZ bendCenter, double radius, double crossRatio, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT_RANGE(startParameter, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT_RANGE(endParameter, knotVector[0], knotVector[knotVector.size() - 1]);
	VALIDATE_ARGUMENT(MathUtils::IsGreaterThan(endParameter, startParameter), "endParameter", "EndParameter must be greater than startParameter.");

	int spanMinIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, startParameter);
	int spanMaxIndex = Polynomials::GetKnotSpanIndex(degree, knotVector, endParameter);

	std::vector<XYZW> updatedControlPoints = controlPoints;
	std::map<int, XYZ> selectedControlPoints;
	for (int i = spanMinIndex; i <= spanMaxIndex - degree - 1; i++)
	{
		XYZ p = const_cast<XYZW&>(updatedControlPoints[i]).ToXYZ(true);
		selectedControlPoints.insert({ i, p });
	}

	for (auto it = selectedControlPoints.begin(); it != selectedControlPoints.end(); ++it)
	{
		int index = it->first;
		XYZ current = it->second;
		XYZ pointOnBendCurve = bendCenter + (current - bendCenter).Normalize() * radius;
		double si = bendCenter.Distance(pointOnBendCurve) / bendCenter.Distance(current);
		double ti = (crossRatio * si) / (1 + (crossRatio - 1) * si);
		XYZ project = bendCenter + (current - bendCenter) * ti;
		updatedControlPoints[index] = XYZW(project, updatedControlPoints[index].GetW());
	}

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = updatedControlPoints;
}

void LNLib::NurbsCurve::ConstraintBasedModification(const LN_NurbsCurve& curve, const std::vector<double>& constraintParams, const std::vector<XYZ>& derivativeConstraints, const std::vector<int>& appliedIndices, const std::vector<int>& appliedDegree, const std::vector<int>& fixedControlPointIndices, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	VALIDATE_ARGUMENT(constraintParams.size() > 0, "constraintParams", "ConstraintParams size must be greater than zero.");
	VALIDATE_ARGUMENT(derivativeConstraints.size() > 0, "derivativeConstraints", "DerivativeConstraints size must be greater than zero.");
	VALIDATE_ARGUMENT(appliedIndices.size() == derivativeConstraints.size(), "appliedIndices", "AppliedIndices size must be equal to derivativeConstraints size.");
	VALIDATE_ARGUMENT(appliedDegree.size() == appliedIndices.size(), "appliedDegree", "AppliedDegree size must be equal to appliedIndices size.");

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

	std::vector<int> remove(size, 1);
	std::vector<int> map(size);

	for (int j = 0; j < size; j++)
	{
		for (int i = 0; i < D.size(); i++)
		{
			if (MathUtils::IsGreaterThan(B[i][j] * B[i][j], 0.0))
			{
				remove[j] = 0;
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
		for (int i = 0; i < B.size(); i++)
		{
			Bopt[i][j] = B[i][map[j]];
		}
	}

	std::vector<std::vector<double>> Bt;
	MathUtils::Transpose(Bopt, Bt);
	std::vector<std::vector<double>> BoptBt = MathUtils::MatrixMultiply(Bopt, Bt);
	std::vector<std::vector<double>> BBt;
	MathUtils::MakeInverse(BoptBt, BBt);

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

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = updatedControlPoints;
}

void LNLib::NurbsCurve::ToClampCurve(const LN_NurbsCurve& curve, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	int m = knotVector.size() - 1;
	double up = knotVector[degree];
	double ump = knotVector[m - degree];

	int upMulti = Polynomials::GetKnotMultiplicity(knotVector, up);
	int umpMulti = Polynomials::GetKnotMultiplicity(knotVector, ump);
	
	LN_NurbsCurve tc;
	int t1 = InsertKnot(curve, up, degree + 1 - upMulti, tc);
	LN_NurbsCurve temp;
	temp = tc;
	int t2 = InsertKnot(temp, ump, degree + 1 - umpMulti, tc);

	std::vector<double> kv = tc.KnotVector;
	std::vector<XYZW> cps = tc.ControlPoints;
	knotVector.resize(kv.size() - t1 - t2, 0.0);
	controlPoints.resize(knotVector.size() - degree - 1);
	for (int i = knotVector.size() - 1; i >= 0; i--)
	{
		knotVector[i] = kv[i + degree];
		if (i < controlPoints.size())
		{
			controlPoints[i] = cps[i + degree];
		}
	}

	result.Degree = tc.Degree;
	result.KnotVector = knotVector;
	result.ControlPoints = controlPoints;

	bool isClosed = IsClosed(curve);
	bool isNowClosed = IsClosed(result);
	if (isClosed && !isNowClosed)
	{
		controlPoints[controlPoints.size() - 1] = controlPoints[0];
		result.ControlPoints = controlPoints;
	}
}

bool LNLib::NurbsCurve::IsPeriodic(const LN_NurbsCurve& curve)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;
	int size = controlPoints.size();

	bool isClamp = IsClamp(curve);
	if (isClamp)
	{
		return false;
	}
	bool isUniform = KnotVectorUtils::IsUniform(knotVector);
	if (!isUniform)
	{
		return false;
	}

	if (size >= degree + degree)
	{
		bool flag = true;
		for (int i = 0; i < degree; i++)
		{
			if (!controlPoints[i].IsAlmostEqualTo(controlPoints[size - degree + i]))
			{
				flag = false;
				break;
			}
		}
		if (flag)
		{
			return true;
		}
	}

	bool isClosed = IsClosed(curve);
	if (!isClosed)
	{
		return false;
	}

	double first = knotVector[0];
	double end = knotVector[knotVector.size() - 1];

	int cFirst = KnotVectorUtils::GetContinuity(curve.Degree, knotVector, first);
	int cEnd = KnotVectorUtils::GetContinuity(curve.Degree, knotVector, end);

	if (cFirst != cEnd)
	{
		return false;
	}

	std::vector<XYZ> fDers = ComputeRationalCurveDerivatives(curve, cFirst, first);
	std::vector<XYZ> eDers = ComputeRationalCurveDerivatives(curve, cEnd, end);

	for (int i = 0; i <= cFirst; i++)
	{
		XYZ currentF = fDers[i];
		XYZ currentE = eDers[i];
		bool hasSameDirection = currentF.Normalize().IsAlmostEqualTo(currentE.Normalize());
		bool hasSameMagnitude = MathUtils::IsAlmostEqualTo(currentF.Length(), currentE.Length());
		if (!hasSameDirection || !hasSameMagnitude)
		{
			return false;
		}
	}
	return true;
}

void LNLib::NurbsCurve::ToUnclampCurve(const LN_NurbsCurve& curve, LN_NurbsCurve& result)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	int n = controlPoints.size() - 1;
	for (int i = 0; i <= degree - 2; i++)
	{
		knotVector[degree - i - 1] = knotVector[degree - i] - (knotVector[n - i + 1] - knotVector[n - i]);
		int k = degree - 1;
		for (int j = i; j >= 0; j--)
		{
			double alpha = (knotVector[degree] - knotVector[k]) / (knotVector[degree + j + 1] - knotVector[k]);
			controlPoints[j] = (controlPoints[j] - alpha * controlPoints[j + 1]) / (1.0 - alpha);
			k = k - 1;
		}
	}
	knotVector[0] = knotVector[1] - (knotVector[n - degree + 2] - knotVector[n - degree + 1]);
	for (int i = 0; i <= degree - 2; i++)
	{
		knotVector[n + i + 2] = knotVector[n + i + 1] + (knotVector[degree + i + 1] - knotVector[degree + i]);
		for (int j = i; j >= 0; j--)
		{
			double alpha = (knotVector[n + 1] - knotVector[n - j]) / (knotVector[n - j + i + 2] - knotVector[n - j]);
			controlPoints[n - j] = (controlPoints[n - j] - (1.0 - alpha) * controlPoints[n - j - 1]) / alpha;
		}
	}
	knotVector[n + degree + 1] = knotVector[n + degree] + (knotVector[2 * degree] - knotVector[2 * degree - 1]);

	result.Degree = degree;
	result.KnotVector = knotVector;
	result.ControlPoints = controlPoints;
}

bool LNLib::NurbsCurve::IsClamp(const LN_NurbsCurve& curve)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	double first = knotVector[0];
	double end = knotVector[knotVector.size() - 1];

	int f = Polynomials::GetKnotMultiplicity(knotVector, first);
	int e = Polynomials::GetKnotMultiplicity(knotVector, end);
	if (f == degree + 1 && e == degree + 1)
	{
		return true;
	}
	return false;
}

bool LNLib::NurbsCurve::IsLinear(const LN_NurbsCurve& curve)
{
	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	int size = controlPoints.size();
	if (size == 2)
		return true;
	XYZ start = controlPoints[0].ToXYZ(true);
	XYZ end =  controlPoints[size - 1].ToXYZ(true);
	if ((end - start).IsAlmostEqualTo(XYZ()))
		return false;
	XYZ dir = (end - start).Normalize();
	for (int i = 0; i < size - 1; i++)
	{
		XYZ current = controlPoints[i].ToXYZ(true);
		XYZ next = controlPoints[i + 1].ToXYZ(true);
		if (!((next - current).Normalize()).IsAlmostEqualTo(dir))
		{
			return false;
		}
	}

	auto map = KnotVectorUtils::GetInternalKnotMultiplicityMap(knotVector);
	for (auto it = map.begin(); it != map.end(); it++)
	{
		double u = it->first;
		XYZ cp = GetPointOnCurve(curve, u);
		XYZ cp2s = (cp - start).Normalize();
		XYZ cp2e = (cp - end).Normalize();

		if (!cp2s.IsAlmostEqualTo(dir))
			return false;
		if (!cp2e.IsAlmostEqualTo(dir.Negative()))
			return false;
	}
	return true;
}

bool LNLib::NurbsCurve::IsArc(const LN_NurbsCurve& curve, LN_ArcInfo& arcInfo)
{
	if (IsLinear(curve))
	{
		return false;
	}

	int degree = curve.Degree;
	std::vector<double> knotVector = curve.KnotVector;
	std::vector<XYZW> controlPoints = curve.ControlPoints;

	double first = knotVector[0];
	double end = knotVector[knotVector.size() - 1];

	XYZ P0 = GetPointOnCurve(curve, first);
	double param = IsClosed(curve)?0.5 * first + 0.5 * end : end;
	XYZ P1 = GetPointOnCurve(curve,0.5 * first + 0.5 * param);
	XYZ P2 = GetPointOnCurve(curve, param);

	XYZ v1 = P1 - P0;
	XYZ v2 = P2 - P0;
	double v1v1, v2v2, v1v2;
	v1v1 = v1.DotProduct(v1);
	v2v2 = v2.DotProduct(v2);
	v1v2 = v1.DotProduct(v2);

	double base = 0.5 / (v1v1 * v2v2 - v1v2 * v1v2);
	double k1 = base * v2v2 * (v1v1 - v1v2);
	double k2 = base * v1v1 * (v2v2 - v1v2);
	XYZ center = P0 + v1 * k1 + v2 * k2;
	double radius = center.Distance(P0);

	std::vector<XYZ> tessellatedPoints;
	std::vector<double> correspondingKnots;
	EquallyTessellate(curve, tessellatedPoints, correspondingKnots);
	for (int i = 0; i < tessellatedPoints.size(); i++)
	{
		XYZ point = tessellatedPoints[i];
		double d = point.Distance(center);
		if (!MathUtils::IsAlmostEqualTo(d, radius))
		{
			return false;
		}
	}
	arcInfo.Radius = radius;
	arcInfo.Center = center;
	return true;
}

double LNLib::NurbsCurve::ApproximateLength(const LN_NurbsCurve& curve, IntegratorType type)
{
	if (IsLinear(curve))
	{
		std::vector<XYZW> controlPoints = curve.ControlPoints;
		XYZ startPoint = controlPoints[0].ToXYZ(true);
		XYZ endPoint = controlPoints[controlPoints.size() - 1].ToXYZ(true);
		return startPoint.Distance(endPoint);
	}

	LN_NurbsCurve reCurve;
	Reparametrize(curve, 0.0, 1.0, reCurve);

	int degree = reCurve.Degree;
	std::vector<double> knotVector = reCurve.KnotVector;
	std::vector<XYZW> controlPoints = reCurve.ControlPoints;

	double length = 0.0;
	switch (type)
	{
		case IntegratorType::Simpson:
		{
			double start = knotVector[0];
			double end = knotVector[knotVector.size() - 1];
			FirstDerivativeLengthFunction function;
			double simpson = Integrator::Simpson(function, (void*)&reCurve, start, end);
			length = CalculateLengthBySimpson(function, reCurve, start, end, simpson, Constants::DistanceEpsilon);
			break;
		}
		case IntegratorType::GaussLegendre:
		{
			// Strongly recommend read this blog:
			// https://raphlinus.github.io/curves/2018/12/28/bezier-arclength.html

			std::vector<LN_NurbsCurve> bezierCurves =  DecomposeToBeziers(reCurve);
			for (int i = 0; i < bezierCurves.size(); i++)
			{
				LN_NurbsCurve bezierCurve = bezierCurves[i];

				std::vector<double> bKnots = bezierCurve.KnotVector;
				double a = bKnots[0];
				double b = bKnots[bKnots.size() - 1];
				double coefficient = (b - a) / 2.0;

				double bLength = 0.0;
				std::vector<double> abscissae = Integrator::GaussLegendreAbscissae;
				int size = abscissae.size();
				for (int i = 0; i < size; i++)
				{
					double t = coefficient * abscissae[i] + (a + b)/2.0;
					double derLength = ComputeRationalCurveDerivatives(bezierCurve, 1, t)[1].Length();
					if (std::isnan(derLength))
						derLength = 0.0;
					bLength += Integrator::GaussLegendreWeights[i] * derLength;
				}
				bLength = coefficient * bLength;
				length += bLength;
			}
			break;
		}
		case IntegratorType::Chebyshev:
		{
			std::vector<double> series = Integrator::ChebyshevSeries();
			for (int i = degree; i < controlPoints.size(); i++) 
			{
				double a = knotVector[i];
				double b = knotVector[i + 1];
				FirstDerivativeLengthFunction function;
				length += Integrator::ClenshawCurtisQuadrature(function, (void*)&reCurve, a, b, series);
			}
			break;
		}
		default:
			 break;
	}
	return length;
}

double LNLib::NurbsCurve::GetParamOnCurve(const LN_NurbsCurve& curve, double givenLength, IntegratorType type)
{
	std::vector<double> knotVector = curve.KnotVector;
	double start = knotVector[0];
	double end = knotVector[knotVector.size() - 1];

	double totalLength = ApproximateLength(curve, type);
	if (MathUtils::IsLessThan(totalLength, givenLength, Constants::DistanceEpsilon))
	{
		return end;
	}

	for (int i = 0; i < knotVector.size(); i++)
	{
		LN_NurbsCurve left;
		LN_NurbsCurve right;
		double knot = knotVector[i];
		bool isSplited = SplitAt(curve, knot, left, right);
		if (!isSplited) continue;
		double length = ApproximateLength(left, type);
		if (MathUtils::IsAlmostEqualTo(length, givenLength))
		{
			return knot;
		}
		if (MathUtils::IsGreaterThan(length, givenLength))
		{
			end = knot;
			break;
		}
	}
	double param = GetParamByLength(curve, start, end, givenLength, type);
	return param;
}

std::vector<double> LNLib::NurbsCurve::GetParamsOnCurve(const LN_NurbsCurve& curve, double givenLength, IntegratorType type)
{
	std::vector<double> result;

	std::vector<double> knotVector = curve.KnotVector;
	double end = knotVector[knotVector.size() - 1];
	
	double param = GetParamOnCurve(curve, givenLength, type);
	while (!MathUtils::IsAlmostEqualTo(param, end))
	{
		result.emplace_back(param);

		LN_NurbsCurve left;
		LN_NurbsCurve right;
		bool isSplited = SplitAt(curve, param, left, right);
		if (!isSplited) continue;
		param = GetParamOnCurve(right, givenLength, type);
	}
	return result;
}

std::vector<LNLib::XYZ> LNLib::NurbsCurve::Tessellate(const LN_NurbsCurve& curve)
{
	int degree = curve.Degree;
	const std::vector<double>& knotVector = curve.KnotVector;
	const std::vector<XYZW>& controlPoints = curve.ControlPoints;

	if (degree == 1)
	{
		return ControlPointsUtils::ToXYZ(controlPoints);
	}

	std::vector<double> unqiueKnotVector = knotVector;
	unqiueKnotVector.erase(std::unique(unqiueKnotVector.begin(), unqiueKnotVector.end()), unqiueKnotVector.end());
	
	std::vector<double> parameters;
	parameters.emplace_back(knotVector[0]);
	for (int i = 0; i < unqiueKnotVector.size() - 1; i++)
	{
		double current = unqiueKnotVector[i];
		double next = unqiueKnotVector[i + 1];
		if (MathUtils::IsLessThanOrEqual(std::abs(next - current), Constants::DoubleEpsilon))
		{
			continue;
		}
		LN_NurbsCurve segment;
		bool isSeg = Segment(curve, current, next, segment);
		if (isSeg)
		{
			bool isLinear = IsLinear(segment);
			if (isLinear)
			{
				parameters.emplace_back(current);
				parameters.emplace_back(next);
				continue;
			}
		}
		TessellateCore(curve, current, next, parameters);
	}
	parameters.emplace_back(knotVector[knotVector.size()-1]);
	std::set<double> st(parameters.begin(), parameters.end());
		
	std::vector<XYZ> points;
	std::set<double>::iterator it = st.begin();
	for (it = st.begin(); it != st.end(); it++)
	{
		XYZ point = GetPointOnCurve(curve, *it);
		points.emplace_back(point);
	}
	return points;
}

