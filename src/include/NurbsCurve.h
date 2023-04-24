#pragma once

#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib
{
	class XYZ;
	class LNLIB_EXPORT NurbsCurve
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page124
		/// Algorithm A4.1
		/// Compute point on rational B-spline curve.
		/// </summary>
		static void GetPointOnCurve(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, const std::vector<double> weights, double paramT, XYZ& point);

		/// <summary>
		/// The NURBS Book 2nd Edition Page127
		/// Algorithm A4.2
		/// Compute C(paramT) derivatives from Cw(paramT) deraivatives.
		/// </summary>
		static void ComputeRationalCurveDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints,  const std::vector<double> weights, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives);
	};
}


