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
