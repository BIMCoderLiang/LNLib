/*
 * Author:
 * 2023/10/28 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#include "ControlPointsUtils.h"
#include "XYZ.h"
#include "XYZW.h"
#include "LNLibExceptions.h"
#include <algorithm>

using namespace LNLib;

std::vector<XYZ> LNLib::ControlPointsUtils::ToXYZ(const std::vector<XYZW>& weightedControlPoints)
{
	std::vector<XYZ> result;
	for (int i = 0; i < weightedControlPoints.size(); i++)
	{
		result.emplace_back(const_cast<XYZW&>(weightedControlPoints[i]).ToXYZ(true));
	}
	return result;
}
