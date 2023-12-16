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

std::vector<std::vector<XYZ>> LNLib::ControlPointsUtils::ToXYZ(const std::vector<std::vector<XYZW>>& points)
{
	int row = points.size();
	int column = points[0].size();

	std::vector<std::vector<XYZ>> result(row, std::vector<XYZ>(column));
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			result[i][j] = const_cast<XYZW&>(points[i][j]).ToXYZ(true);
		}
	}
	return result;
}

std::vector<std::vector<XYZW>> LNLib::ControlPointsUtils::ToXYZW(const std::vector<std::vector<XYZ>>& points)
{
	int row = points.size();
	int column = points[0].size();

	std::vector<std::vector<XYZW>> result(row, std::vector<XYZW>(column));;
	for (int i = 0; i < row; i++)
	{
		for (int j = 0; j < column; j++)
		{
			result[i][j] = XYZW(const_cast<XYZ&>(points[i][j]), 1);
		}
	}
	return result;
}

std::vector<std::vector<XYZW>> LNLib::ControlPointsUtils::Multiply(const std::vector<std::vector<XYZW>>& points, const std::vector<std::vector<double>>& coefficient)
{
	std::vector<std::vector<XYZW>> result(coefficient.size(), std::vector<XYZW>(points[0].size()));
	for (int i = 0; i < coefficient.size(); i++)
	{
		for (int j = 0; j < points[0].size(); j++)
		{
			for (int k = 0; k < coefficient[0].size(); k++)
			{
				result[i][j] += coefficient[i][k] * points[k][j];
			}
		}
	}
	return result;
}
