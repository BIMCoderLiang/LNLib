/*
 * Author: 
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "BezierCurve.h"
#include "ValidationUtils.h"
#include "XYZ.h"
#include "XYZW.h"
#include "Matrix4d.h"
#include "Polynomials.h"
#include "Intersection.h"
#include "Interpolation.h"
#include "ValidationUtils.h"
#include "LNLibExceptions.h"
#include <vector>

using namespace LNLib;

XYZ LNLib::BezierCurve::GetPointOnQuadraticArc(const XYZW& startPoint, const XYZW& middlePoint, const XYZW& endPoint, double paramT)
{
	VALIDATE_ARGUMENT_RANGE(paramT, 0, 1);

	double w0 = startPoint.GetW();
	XYZ P0 = XYZ(startPoint.GetWX()/w0, startPoint.GetWY()/w0, startPoint.GetWZ()/w0);
	double w1 = middlePoint.GetW();
	XYZ P1 = XYZ(middlePoint.GetWX() / w0, middlePoint.GetWY() / w0, middlePoint.GetWZ() / w0);
	double w2 = endPoint.GetW();
	XYZ P2 = XYZ(endPoint.GetWX() / w0, endPoint.GetWY() / w0, endPoint.GetWZ() / w0);

	return  ((1 - paramT) * (1 - paramT) * w0 * P0 + 2 * paramT * (1 - paramT) * w1 * P1 + paramT * paramT * w2 * P2) / ((1 - paramT) * (1 - paramT) * w0 + 2 * paramT * (1 - paramT) * w1 + paramT * paramT * w2);
}

bool LNLib::BezierCurve::ComputerMiddleControlPointsOnQuadraticCurve(const XYZ& startPoint, const XYZ& startTangent, const XYZ& endPoint, const XYZ& endTangent, std::vector<XYZW>& controlPoints)
{
	XYZ chord = (endPoint - startPoint).Normalize();
	XYZ nST = const_cast<XYZ&>(startTangent).Normalize();
	XYZ nET = const_cast<XYZ&>(endTangent).Normalize();

	if (nST.IsAlmostEqualTo(chord) ||
		nST.IsAlmostEqualTo(chord.Negative()))
	{
		if (nST.IsAlmostEqualTo(nET) ||
			nST.IsAlmostEqualTo(nET.Negative()))
		{
			controlPoints.emplace_back((startPoint + endPoint) / 2, 1);
			return true;
		}
	}

	double alf1 = 0.0;
	double alf2 = 0.0;

	XYZ R(0, 0, 0);
	double weight = 0.0;
	auto type = Intersection::ComputeRays(startPoint, nST, endPoint, nET, alf1, alf2, R);
	if (type == CurveCurveIntersectionType::Intersecting)
	{
		if (startPoint.IsAlmostEqualTo(R) || endPoint.IsAlmostEqualTo(R))
		{
			R = 0.5 * (startPoint + endPoint);
			bool result = Interpolation::ComputerWeightForRationalQuadraticInterpolation(startPoint, R, endPoint, weight);
			if (!result) return false;
			controlPoints.emplace_back(XYZW(R, weight));
			return true;
		}

		if (MathUtils::IsGreaterThan(alf1, 0.0) && MathUtils::IsLessThan(alf2, 0.0))
		{
			bool result = Interpolation::ComputerWeightForRationalQuadraticInterpolation(startPoint, R, endPoint, weight);
			if (!result) return false;
			controlPoints.emplace_back(XYZW(R, weight));
			return true;
		}
	}
	
	XYZ SE = endPoint - startPoint;
	if (SE.Normalize().IsAlmostEqualTo(startTangent) && SE.Normalize().IsAlmostEqualTo(endTangent))
	{
		R = 0.5 * (startPoint + endPoint);
		bool result = Interpolation::ComputerWeightForRationalQuadraticInterpolation(startPoint, R, endPoint, weight);
		if (!result) return false;
		controlPoints.emplace_back(XYZW(R, weight));
		return true;
	}

	double gamma1 = 0.0;
	double gamma2 = 0.0;

	if (type != CurveCurveIntersectionType::Intersecting)
	{
		gamma1 = 0.5 * SE.Length();
		gamma2 = gamma1;
	}
	else
	{
		XYZ SR = R - startPoint;
		XYZ ER = R - endPoint;

		double theta0 = SR.DotProduct(SE) / (SR.Length() * SE.Length());
		double theta1 = ER.DotProduct(SE) / (ER.Length() * SE.Length());

		double alpha = 2.0 / 3.0;
		gamma1 = 0.5 * SE.Length() / (1.0 + alpha * theta1 + (1 - alpha) * theta0);
		gamma2 = 0.5 * SE.Length() / (1.0 + alpha * theta0 + (1 - alpha) * theta1);
	}

	XYZ R1 = startPoint + gamma1 * startTangent;
	XYZ R2 = endPoint - gamma2 * endTangent;
	XYZ Qk = (gamma1 * R2 + gamma2 * R1) / (gamma1 + gamma2);

	double weight1 = 0.0;
	bool result = Interpolation::ComputerWeightForRationalQuadraticInterpolation(startPoint, R1, Qk, weight1);
	if (!result) return false;

	double weight2 = 0.0;
	result = Interpolation::ComputerWeightForRationalQuadraticInterpolation(Qk, R2, endPoint, weight2);
	if (!result) return false;

	controlPoints.emplace_back(XYZW(R1, weight1));
	controlPoints.emplace_back(XYZW(Qk, 1.0));
	controlPoints.emplace_back(XYZW(R2, weight2));

	return true;
}
