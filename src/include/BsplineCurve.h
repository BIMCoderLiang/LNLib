#pragma once

#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib
{
	class XYZ;
	class XYZW;
	class LNLIB_EXPORT BsplineCurve
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page82
		/// Algorithm A3.1
		/// Compute Bspline curve point
		/// </summary>
		static void GetPointOnCurve(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT, const std::vector<double>& knotVector, XYZ& point);

		/// <summary>
		/// The NURBS Book 2nd Edition Page93
		/// Algorithm A3.2
		/// Compute curve derivatives.
		/// </summary>
		static void ComputeDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives);

		/// <summary>
		/// Compute curve derivatives for XYZW using Algorithm A3.2
		/// </summary>
		static void ComputeDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZW>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZW>& derivatives);

		/// <summary>
		/// The NURBS Book 2nd Edition Page98
		/// Algorithm A3.3
		/// Compute control points of curve derivatives.
		/// </summary>
		static void ComputeControlPointsOfDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, unsigned int derivative, unsigned int min, unsigned int max, std::vector<std::vector<XYZ>>& controlPointsOfDerivative);

		/// <summary>
		/// The NURBS Book 2nd Edition Page99
		/// Algorithm A3.4
		/// Compute curve detivatives.
		/// </summary>
		static void ComputeDerivativesByAllBasisFunctions(unsigned int degree, const std::vector<double>& knotVector, const std::vector<XYZ>& controlPoints, double paramT, unsigned int derivative, std::vector<XYZ>& derivatives);
	};
}


