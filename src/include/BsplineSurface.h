#pragma once


#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib
{
	class XYZ;
	class XYZW;
	class LNLIB_EXPORT BsplineSurface
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page103
		/// Algorithm A3.5
		/// Compute surface point.
		/// </summary>
		static void GetPointOnSurface(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, XYZ& point);

		/// <summary>
		/// The NURBS Book 2nd Edition Page111
		/// Algorithm A3.6
		/// Compute surface derivatives.
		/// </summary>
		static void ComputeDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives);

		/// <summary>
		/// Compute surface derivatives for XYZW using A3.6.
		/// </summary>
		static void ComputeDerivatives(const std::vector<std::vector<XYZW>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, std::vector<std::vector<XYZW>>& derivatives);

		/// <summary>
		/// The NURBS Book 2nd Edition Page114.
		/// Algorithm A3.7
		/// Compute control points of derivative surfaces.
		/// </summary>
		static void ComputeControlPointsOfDerivatives(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, unsigned int minU, unsigned int maxU, unsigned int minV, unsigned int maxV, std::vector<std::vector<std::vector<std::vector<XYZ>>>> & controlPointsOfDerivative);

		/// <summary>
		/// The NURBS Book 2nd Edition Page115.
		/// Algorithm A3.8
		/// Compute surface derivatives.
		/// </summary>
		static void ComputeDerivativesByAllBasisFunctions(const std::vector<std::vector<XYZ>>& controlPoints, const std::vector<double>& knotVectorU, const std::vector<double>& knotVectorV, unsigned int degreeU, unsigned int degreeV, double paramU, double paramV, unsigned int derivative, std::vector<std::vector<XYZ>>& derivatives);
	};
}


