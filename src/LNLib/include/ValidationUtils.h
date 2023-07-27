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
#include "MathUtils.h"
#include <vector>

namespace LNLib
{
	class XYZ;
	class XYZW;
	class LNLIB_EXPORT ValidationUtils
	{
	public:

		static bool IsInRange(double input, double min, double max);

		static bool IsValidBezier(unsigned int degree, unsigned int controlPointsCount);

		static bool IsValidBspline(unsigned int degree, unsigned int knotVectorCount, unsigned int controlPointsCount);

		static bool IsValidNurbs(unsigned int degree, unsigned int knotVectorCount, unsigned int controlPointsCount, unsigned int weightsCount);

		static bool IsValidDegreeReduction(unsigned int degree);

		/// <summary>
		/// The NURBS Book 2nd Edition Page50
		/// Knot Vector is a nondecreasing sequence of real numbers.
		/// </summary>
		static bool IsValidKnotVector(const std::vector<double>& knotVector);

		/// <summary>
		/// The NURBS Book 2nd Edition Page185
		/// TOL = dWmin / (1+abs(Pmax))
		/// </summary>
		static double ComputeCurveModifyTolerance(const std::vector<XYZW>& controlPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page221
		/// Implements Bezier degree reduction and computation of the maximum error.
		/// </summary>
		static double ComputeMaxErrorOfBezierReduction(unsigned int degree, const std::vector<XYZW>& currentControlPoints, const std::vector<XYZW>& reductedControlPoints);

		template <typename T>
		static bool IsClosed(const std::vector<T>& controlPoints)
		{
			T first = controlPoints[0];
			T last = controlPoints[controlPoints.size() - 1];
			return first.IsAlmostEqualTo(last);
		}

		///  [0][0]  [0][1] ... ...  [0][m]     ------- v direction
		///  [1][0]  [1][1] ... ...  [1][m]    |
		///    .                               |
		///    .                               u direction
		///    .							   
		///  [n][0]  [n][1] ... ...  [n][m]      
		template <typename T>
		static bool IsClosedU(const std::vector<std::vector<T>>& controlPoints)
		{
			std::vector<std::vector<T>> transposed;
			MathUtils::Transpose(controlPoints, transposed);

			for (int i = 0; i < transposed.size(); i++)
			{
				std::vector<T> row = transposed[i];
				bool rowResult = IsClosed(row);
				if (!rowResult)
				{
					return false;
				}
			}
			return true;
		}

		template <typename T>
		static bool IsClosedV(const std::vector<std::vector<T>>& controlPoints)
		{
			for (int i = 0; i < controlPoints.size(); i++)
			{
				std::vector<T> row = controlPoints[i];
				bool rowResult = IsClosed(row);
				if (!rowResult)
				{
					return false;
				}
			}
			return true;
		}
	};
}

