#include "BezierSurface.h"
#include "BezierCurve.h"
#include "XYZ.h"
#include "UV.h"

using namespace LNLib;

void BezierSurface::GetPointOnSurfaceByDeCasteljau(const std::vector<std::vector<XYZ>>& controlPoints, unsigned int degreeU, unsigned int degreeV, UV uv, XYZ& point)
{
	
	if (degreeU <= degreeV)
	{		
		std::vector<XYZ> temp;
		temp.resize(degreeV + 1);

		for (unsigned int j = 0; j <= degreeV; j++)
		{
			BezierCurve::GetPointOnCurveByDeCasteljau(controlPoints[j], degreeU, uv.GetU(), temp[j]);
		}
		BezierCurve::GetPointOnCurveByDeCasteljau(temp, degreeV, uv.GetV(), point);
	}
	else
	{
		std::vector<XYZ> temp;
		temp.resize(degreeU + 1);

		int size = static_cast<int>(controlPoints.size());

		for (unsigned int i = 0; i <= degreeU; i++)
		{
			std::vector<XYZ> column;
			column.resize(size);

			for (int k = 0; k < size; k++)
			{
				column[k] = (controlPoints[k][i]);
			}
			BezierCurve::GetPointOnCurveByDeCasteljau(column, degreeV, uv.GetV(), temp[i]);
		}
		BezierCurve::GetPointOnCurveByDeCasteljau(temp, degreeU, uv.GetU(), point);
	}
}
