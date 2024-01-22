/*
 * Author:
 * 2024/01/08 - Yuqing Liang (BIMCoder Liang)
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
	class LNLIB_EXPORT Integrator
	{

	public:
		static double Simpson(double start, double end, double startTangentLength, double middleTangentLength, double endTangentLength);

		static double Simpson(double start, double end, std::vector<double> odds, std::vector<double> evens, double delta);

		/// <summary>
		/// According to https://github.com/Pomax/bezierjs
		/// Order n = 24.
		/// </summary>
		static const std::vector<double> GaussLegendreAbscissae;
		static const std::vector<double> GaussLegendreWeights;

		static std::vector<double> ClenshawCurtisQuadratureWeights(int degree);
	};

}


