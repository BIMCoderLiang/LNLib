#pragma once

#include "LNLibDefinitions.h"
#include "Polynomials.h"
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
		/// Compute Bspline curve point.
		/// </summary>
		static void GetPointOnCurve(const std::vector<XYZ>& controlPoints, unsigned int degree, double paramT, const std::vector<double>& knotVector, XYZ& point);

		/// <summary>
		/// The NURBS Book 2nd Edition Page88
		/// Compute the continuity.
		/// </summary>
		static int GetContinuity(unsigned int degree, const std::vector<double>& knotVector, double knot);

		/// <summary>
		/// The NURBS Book 2nd Edition Page93
		/// Algorithm A3.2
		/// Compute curve derivatives. (Usually Use)
		/// </summary>
		template<typename T>
		static void ComputeDerivatives(unsigned int degree, const std::vector<double>& knotVector, const std::vector<T>& controlPoints, double paramT, unsigned int derivative, std::vector<T>& derivatives)
		{
			derivatives.resize(derivative + 1);

			int du = std::min(derivative, degree);
			int n = static_cast<int>(controlPoints.size() - 1);
			int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);

			std::vector<std::vector<double>> nders;
			Polynomials::BasisFunctionsDerivatives(spanIndex, degree, paramT, du, knotVector, nders);

			for (int k = 0; k <= du; k++)
			{
				for (unsigned int j = 0; j <= degree; j++)
				{
					derivatives[k] += nders[k][j] * controlPoints[spanIndex - degree + j];
				}
			}
		}

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
		template<typename T>
		static void ComputeDerivativesByAllBasisFunctions(unsigned int degree, const std::vector<double>& knotVector, const std::vector<T>& controlPoints, double paramT, unsigned int derivative, std::vector<T>& derivatives)
		{
			derivatives.resize(derivative + 1);

			int du = std::min(derivative, degree);
			int n = static_cast<int>(controlPoints.size() - 1);
			int spanIndex = Polynomials::GetKnotSpanIndex(n, degree, paramT, knotVector);

			std::vector<std::vector<double>> allBasisFunctions;
			allBasisFunctions.resize(degree + 1);
			for (int i = 0; i <= static_cast<int>(degree); i++)
			{
				allBasisFunctions[i].resize(degree + 1);
			}
			for (int i = 0; i <= static_cast<int>(degree); i++)
			{
				std::vector<double> basisFunctions;
				Polynomials::BasisFunctions(spanIndex, i, paramT, knotVector, basisFunctions);
				allBasisFunctions[i] = basisFunctions;
			}

			std::vector<std::vector<T>> controlPointsOfDerivative;
			ComputeControlPointsOfDerivatives(degree, knotVector, controlPoints, du, spanIndex - degree, spanIndex, controlPointsOfDerivative);

			for (int k = 0; k <= du; k++)
			{
				for (int j = 0; j <= static_cast<int>(degree) - k; j++)
				{
					derivatives[k] += allBasisFunctions[j][degree - k] * controlPointsOfDerivative[k][j];
				}
			}
		}
	};
}


