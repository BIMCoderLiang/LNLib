/*
 * Author:
 * 2023/07/04 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
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
		static double GetTotalChordLength(const std::vector<XYZ>& throughPoints);

		static std::vector<double> GetChordParameterization(const std::vector<XYZ>& throughPoints);

		static void ComputeKnotVector(unsigned int degree, const std::vector<XYZ>& throughPoints, std::vector<double>& knotVector);

		static bool LUDecomposition(const std::vector<std::vector<double>>& matrix, std::vector<std::vector<double>>& matrixL, std::vector<std::vector<double>>& matrixU);

		static std::vector<double> ForwardSubstitution(const std::vector<std::vector<double>>& matrixL, const std::vector<double>& column);

		static std::vector<double> BackwardSubstitution(const std::vector<std::vector<double>>& matrixU, const std::vector<double>& column);
	};
}