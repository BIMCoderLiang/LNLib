/*
 * Author:
 * 2023/10/28 - Yuqing Liang (BIMCoder Liang)
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
	class XYZ;
	class XYZW;
	class LNLIB_EXPORT ControlPointsUtils
	{
	public:
		static std::vector<XYZ> ToXYZ(const std::vector<XYZW>& weightedControlPoints);
		
		static std::vector<std::vector<XYZ>> ToXYZ(const std::vector<std::vector<XYZW>>& points);

		static std::vector<std::vector<XYZW>> ToXYZW(const std::vector<std::vector<XYZ>>& points);

		static std::vector<std::vector<XYZW>> Multiply(const std::vector<std::vector<XYZW>>& points, const std::vector<std::vector<double>>& coefficient);
	};

}



