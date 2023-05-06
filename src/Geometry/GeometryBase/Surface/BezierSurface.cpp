#include "BezierSurface.h"
#include "BezierCurve.h"
#include "XYZ.h"

using namespace LNLib;

void BezierSurface::GetPointOnSurfaceByDeCasteljau(const std::vector<std::vector<XYZ>>& controlPoints, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, XYZ& point)
{
	
	if (degreeU <= degreeV)
	{		
		std::vector<XYZ> temp;
		temp.resize(degreeV + 1);

		for (unsigned int j = 0; j <= degreeV; j++)
		{
			BezierCurve::GetPointOnCurveByDeCasteljau(controlPoints[j], degreeU, paramU, temp[j]);
		}
		BezierCurve::GetPointOnCurveByDeCasteljau(temp, degreeV, paramV, point);
	}
	else
	{
		std::vector<XYZ> temp;
		temp.resize(degreeU + 1);

		for (unsigned int i = 0; i <= degreeU; i++)
		{
			std::vector<XYZ> column;
			column.resize(controlPoints.size());

			for (unsigned int k = 0; k < controlPoints.size(); k++)
			{
				column[k] = (controlPoints[k][i]);
			}
			BezierCurve::GetPointOnCurveByDeCasteljau(column, degreeV, paramV, temp[i]);
		}
		BezierCurve::GetPointOnCurveByDeCasteljau(temp, degreeU, paramU, point);
	}
}
