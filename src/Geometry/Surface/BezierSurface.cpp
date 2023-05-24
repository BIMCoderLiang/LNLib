#include "BezierSurface.h"
#include "BezierCurve.h"
#include "XYZ.h"
#include "UV.h"

using namespace LNLib;

void BezierSurface::GetPointOnSurfaceByDeCasteljau(const std::vector<std::vector<XYZ>>& controlPoints, unsigned int n, unsigned int m, UV uv, XYZ& point)
{
	if (n <= m)
	{		
		std::vector<XYZ> temp;
		temp.resize(m + 1);

		for (unsigned int j = 0; j <= m; j++)
		{
			BezierCurve::GetPointOnCurveByDeCasteljau(controlPoints[j], n, uv.GetU(), temp[j]);
		}
		BezierCurve::GetPointOnCurveByDeCasteljau(temp, m, uv.GetV(), point);
	}
	else
	{
		std::vector<XYZ> temp;
		temp.resize(n + 1);

		int size = static_cast<int>(controlPoints.size());

		for (unsigned int i = 0; i <= n; i++)
		{
			std::vector<XYZ> column;
			column.resize(size);

			for (int k = 0; k < size; k++)
			{
				column[k] = (controlPoints[k][i]);
			}
			BezierCurve::GetPointOnCurveByDeCasteljau(column, m, uv.GetV(), temp[i]);
		}
		BezierCurve::GetPointOnCurveByDeCasteljau(temp, n, uv.GetU(), point);
	}
}
