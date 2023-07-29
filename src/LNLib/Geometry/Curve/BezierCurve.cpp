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
#include "Polynomials.h"
#include "ValidationUtils.h"
#include <vector>

using namespace LNLib;

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

XYZW LNLib::BezierCurve::ComputerMiddleControlPointOnQuadraticCurve(const XYZ& startPoint, const XYZ& startTangent, const XYZ& endPoint, const XYZ& endTangent)
{
	XYZ middlePoint = XYZ(0, 0, 0);
	double mw = 1.0;

	// to be continued...
	return XYZW(middlePoint, mw);
}
