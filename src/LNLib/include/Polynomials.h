/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once

#include "LNLibDefinitions.h"
#include <vector>
#include <unordered_map>

namespace LNLib
{
	class UV;
	class LNLIB_EXPORT Polynomials
	{

	public:
		
		/// <summary>
		/// The NURBS Book 2nd Edition Page20
		/// Algorithm A1.1
		/// Compute point on power basis curve.
		/// </summary>
		static double Horner(const std::vector<double>& coefficients, int degree, double paramT);

		/// <summary>
		/// The NURBS Book 2nd Edition Page7
		/// Algorithm A1.2
		/// compute the value of a Berstein polynomial.
		/// </summary>
		static double Bernstein(int index, int degree, double paramT);

		/// <summary>
		/// The NURBS Book 2nd Edition Page21
		/// Algorithm A1.3
		/// Compute All nth-degree Berstein polynomials.
		/// </summary>
		static std::vector<double> AllBernstein(int degree, double paramT);

		/// <summary>
		/// The NURBS Book 2nd Edition Page36
		/// Algorithm A1.6
		/// Compute point on a power basis surface.
		/// coefficients with (n+1) * (m+1)
		///  
		///  [0][0]  [0][1] ... ...  [0][m]     ------- v direction
		///  [1][0]  [1][1] ... ...  [1][m]    |
		///    .                               |
		///    .                               u direction
		///    .
		///  [n][0]  [n][1] ... ...  [n][m]      
		/// 
		/// </summary>
		static double Horner(const std::vector<std::vector<double>>& coefficients, unsigned int n, unsigned int m, UV& uv);

		/// <summary>
		/// The NURBS Book 2nd Edition Page68
		/// Algorithm A2.1
		/// Get the knot span index.
		/// </summary>
		static int GetKnotSpanIndex(unsigned int n, unsigned int degree, double paramT, const std::vector<double>& knotVector);

		/// <summary>
		/// The NURBS Book 2nd Edition Page152
		/// Get the knot multiplicity.
		/// </summary>
		static int GetKnotMultiplicity(double knot, const std::vector<double>& knotVector);

		static void GetInsertedKnotElement(const std::vector<double> knotVector0, const std::vector<double> knotVector1, std::vector<double>& insertElements0, std::vector<double>& insertElements1);

		/// <summary>
		/// The NURBS Book 2nd Edition Page70
		/// Algorithm A2.2
		/// Compute the nonvanishing basis functions.
		/// </summary>
		static void BasisFunctions(unsigned int spanIndex, unsigned int degree, double paramT, const std::vector<double>& knotVector, std::vector<double>& basisFunctions);

		/// <summary>
		/// The NURBS Book 2nd Edition Page72
		/// Algorithm A2.3
		/// Compute nonzero basis functions and their derivative.
		/// </summary>
		static void BasisFunctionsDerivatives(unsigned int spanIndex, unsigned int degree, double paramT, unsigned int derivative, const std::vector<double>& knotVector, std::vector<std::vector<double>>& derivatives);

		/// <summary>
		/// The NURBS Book 2nd Edition Page74
		/// Algorithm A2.4
		/// Compute a single basis function.
		/// </summary>
		static double OneBasisFunction(unsigned int index, unsigned int degree, const std::vector<double>& knotVector, double paramT);

		/// <summary>
		/// The NURBS Book 2nd Edition Page76
		/// Algorithm A2.5
		/// Compute a single basis function and its derivative.
		/// </summary>
		static void OneBasisFunctionDerivative(unsigned int index, unsigned int degree, const std::vector<double>& knotVector, double paramT, unsigned int derivative, std::vector<double>& derivatives);

		/// <summary>
		/// The NURBS Book 2nd Edition Page269
		/// Algorithm A6.1
		/// Compute pth degree Bezier matrix.
		/// </summary>
		static void BezierToPowerMatrix(unsigned int degree, std::vector<std::vector<double>>& matrix);

		/// <summary>
		/// The NURBS Book 2nd Edition Page275
		/// Algorithm A6.2
		/// Compute inverse of pth-degree Bezier matrix.
		/// </summary>
		static void PowerToBezierMatrix(unsigned int degree, const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& inverseMatrix);
	};

}



