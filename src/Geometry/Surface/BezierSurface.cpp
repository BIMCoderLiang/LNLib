/**
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#include "BezierSurface.h"
#include "BezierCurve.h"
#include "XYZ.h"
#include "UV.h"

using namespace LNLib;

void BezierSurface::GetPointOnSurfaceByDeCasteljau(const std::vector<std::vector<XYZ>>& controlPoints, unsigned int n, unsigned int m, UV uv, XYZ& point)
{
	std::vector<XYZ> temp;
	temp.resize(n + 1);

	for (unsigned int i = 0; i <= n; i++)
	{
		BezierCurve::GetPointOnCurveByDeCasteljau(controlPoints[i], m, uv.GetV(), temp[i]);
	}
	BezierCurve::GetPointOnCurveByDeCasteljau(temp, n, uv.GetU(), point);
}
