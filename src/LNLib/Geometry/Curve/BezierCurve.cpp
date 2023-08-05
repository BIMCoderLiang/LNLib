/*
 * Author: 
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "BezierCurve.h"
#include "ValidationUtils.h"
#include "XYZ.h"
#include "XYZW.h"
#include "Matrix4d.h"
#include "Polynomials.h"
#include "Intersection.h"
#include "ValidationUtils.h"
#include <vector>

using namespace LNLib;

namespace LNLib
{
	bool ComputerWeightForRationalQuadraticInterpolation(const XYZ& startPoint, const XYZ& endPoint, const XYZ& middleControlPoint, double weight)
	{
		XYZ R = middleControlPoint;
		XYZ normal = (R - startPoint).Normalize().CrossProduct((R - endPoint).Normalize());
		XYZ M = (startPoint + endPoint) / 2;

		double angle1 = (M - startPoint).Normalize().AngleTo((R - startPoint).Normalize());
		double halfangle1 = angle1 / 2;
		Matrix4d matrix = Matrix4d::CreateRotation(normal, halfangle1);
		XYZ t1 = matrix.OfVector((M - startPoint).Normalize());

		double a1 = 0.0, a2 = 0.0;
		XYZ S1(0, 0, 0);
		CurveCurveIntersectionType type1 = Intersection::ComputeRays(startPoint, t1, M, (R - M).Normalize(), a1, a2, S1);

		double angle2 = (M - endPoint).Normalize().AngleTo((R - endPoint).Normalize());
		double halfangle2 = angle2 / 2;
		matrix = Matrix4d::CreateRotation(normal, halfangle2);
		XYZ t2 = matrix.OfVector((M - endPoint).Normalize());

		XYZ S2(0, 0, 0);
		CurveCurveIntersectionType type2 = Intersection::ComputeRays(endPoint, t2, M, (R - M).Normalize(), a1, a2, S2);

		if (type1 == CurveCurveIntersectionType::Intersecting &&
			type2 == CurveCurveIntersectionType::Intersecting)
		{
			XYZ S = (S1 + S2) / 2;
			double s = (S - M).Length() / (R - M).Length();
			if ((S - M).Normalize().IsAlmostEqualTo((R - M).Normalize().Negative()))
			{
				s = -s;
			}
			weight = s / (1 - s);
			return true;
		}
		return false;
	}
}

XYZ BezierCurve::GetPointOnCurveByBernstein(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT)
{
	std::vector<double> bernsteinArray;
	Polynomials::AllBernstein(degree, paramT, bernsteinArray);

	XYZ temp(0,0,0);
	for (unsigned int k = 0; k <= degree; k++)
	{
		temp += bernsteinArray[k] * controlPoints[k];
	}
	return temp;
}

XYZ BezierCurve::GetPointOnCurveByDeCasteljau(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT)
{
	std::vector<XYZ> temp;
	temp.resize(degree + 1);

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
	return temp[0];
}

XYZ LNLib::BezierCurve::GetPointOnQuadraticArc(const XYZW& startPoint, const XYZW& middlePoint, const XYZW& endPoint, double paramT)
{
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

	double a1 = 0.0, a2 = 0.0;
	XYZ R(0, 0, 0);
	CurveCurveIntersectionType type = Intersection::ComputeRays(startPoint, nST, endPoint, nET, a1, a2, R);
	if (type == CurveCurveIntersectionType::Intersecting)
	{
		double d1 = startPoint.Distance(R);
		double d2 = endPoint.Distance(R);
		if (MathUtils::IsAlmostEqualTo(d1, d2))
		{
			double w = ((startPoint + endPoint) / 2).Distance(endPoint) / (R.Distance(endPoint));
			controlPoints.emplace_back(R, w);
			return true;
		}
		else
		{
			double gamma1 = (R - startPoint)[0] / startTangent[0];
			double gamma2 = (R - endPoint)[0] / endTangent[0];
			if (MathUtils::IsGreaterThan(gamma1, 0.0) && MathUtils::IsLessThan(gamma2, 0.0))
			{
				double w = 0.0;
				if (ComputerWeightForRationalQuadraticInterpolation(startPoint, R, endPoint, w))
				{
					controlPoints.emplace_back(R, w);
					return true;
				}
				else
				{
					return false;
				}
			}
			else
			{
				if (nST.IsAlmostEqualTo(nET) ||
					nST.IsAlmostEqualTo(nET.Negative()))
				{
					gamma1 = gamma2 = 0.5 * startPoint.Distance(endPoint);
				}
				else
				{
					double thetak1 = (endPoint - startPoint).Normalize().AngleTo(nST);
					if (MathUtils::IsGreaterThan(thetak1, Constants::Pi / 2))
					{
						thetak1 = Constants::Pi - thetak1;
					}
					double thetak = (endPoint - startPoint).Normalize().AngleTo(nET);
					if (MathUtils::IsGreaterThan(thetak, Constants::Pi / 2))
					{
						thetak = Constants::Pi - thetak;
					}
					gamma1 = 0.5 * (endPoint.Distance(startPoint)) / (1 + (2 / 3) * cos(thetak) + (1 - (2 / 3)) * cos(thetak1));
					gamma2 = 0.5 * (endPoint.Distance(startPoint)) / (1 + (2 / 3) * cos(thetak1) + (1 - (2 / 3)) * cos(thetak));
				}

				XYZ Rkt = startPoint + gamma1 * startTangent;
				XYZ Rkt1 = endPoint - gamma2 * endTangent;
				XYZ Qkt = (gamma1 * Rkt1 + gamma2 * Rkt) / (gamma1 + gamma2);

				double w1 = 0.0;
				bool b1 = ComputerWeightForRationalQuadraticInterpolation(startPoint, Rkt, Qkt, w1);
				double w2 = 0.0;
				bool b2 = ComputerWeightForRationalQuadraticInterpolation(Qkt, Rkt1, endPoint, w2);
				if (b1 && b2)
				{
					controlPoints.emplace_back(Rkt, w1);
					controlPoints.emplace_back(Rkt1, w2);
					return true;
				}
			}
		}
	}
	return false;
}
