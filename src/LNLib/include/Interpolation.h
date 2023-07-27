/*
 * Author:
 * 2023/07/04 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib
{
	class XYZ;
	class LNLIB_EXPORT Interpolation
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page364
		/// The total chord length.
		/// </summary>
		static double GetTotalChordLength(const std::vector<XYZ>& throughPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page365
		/// The chord length parameterization.
		/// </summary>
		static std::vector<double> GetChordParameterization(const std::vector<XYZ>& throughPoints);

		/// <summary>
		/// The NURBS Book 2nd Edition Page365
		/// Technique of averaging.
		/// </summary>
		static void ComputeKnotVector(unsigned int degree, int pointsCount, const std::vector<double> params, std::vector<double>& knotVector);

		/// <summary>
		/// The NURBS Book 2nd Edition Page412
		/// Computes a knot vector ensuring that every knot span has at least one.
		/// </summary>
		static void ComputeKnotVector(unsigned int degree, int pointsCount, int controlPointsCount, const std::vector<double> params, std::vector<double>& knotVector);

		static std::vector<std::vector<double>> MakeInterpolationMatrix(unsigned int degree, int dataCount, const std::vector<double>& params, const std::vector<double>& knotVector);

		static std::vector<XYZ> ComputerMatrixMultiplyPoints(std::vector<std::vector<double>> matrix, std::vector<XYZ> points);

		static std::vector<XYZ> ComputerControlPointsByLUDecomposition(const std::vector<std::vector<double>>& matrix, const std::vector<XYZ>& data);

		static void ComputerKnotVectorForTangents(unsigned int degree, const std::vector<double>& params, const std::vector<int>& derivativeIndices, std::vector<double>& knotVector);

		/// <summary>
		/// The NURBS Book 2nd Edition Page377
		/// Algorithm A9.3
		/// Compute paramters for global surface interpolation.
		/// </summary>
		static void GetSurfaceMeshParameterization(const std::vector<std::vector<XYZ>>& throughPoints, std::vector<double>& paramVectorU, std::vector<double>& paramVectorV);

		/// <summary>
		/// The NURBS Book 2nd Edition Page386
		/// Computer tangent of each through point (at least five points).
		/// </summary>
		static bool ComputerTangent(const std::vector<XYZ>& throughPoints, std::vector<XYZ>& tangents);
	};
}