#pragma once

#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib
{
	class XYZ;
	class LNLIB_EXPORT BezierCurve
	{

	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page22
		/// Algorithm A1.4
		/// Compute point on Bezier curve.
		/// </summary>
		static void GetPointOnCurveByBernstein(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT, XYZ& point);

		/// <summary>
		/// The NURBS Book 2nd Edition Page24
		/// Algorithm A1.5
		/// Compute point on Bezier curve using deCasteljau.
		/// </summary>
		static void GetPointOnCurveByDeCasteljau(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT, XYZ& point);

	};
}


