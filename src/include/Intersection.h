/*
 * Author:
 * 2023/06/23 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license license that can be found in
 * the LICENSE file.
 */

#pragma once

#include "LNLibDefinitions.h"

namespace LNLib
{
	enum class CurveCurveIntersectionType : int
	{
		Intersecting = 0,
		Parallel = 1,
		Coincident = 2,
		Skew = 3,
	};

	class XYZ;
	class LNLIB_EXPORT Intersection
	{

	public:
		static CurveCurveIntersectionType ComputeRays(const XYZ& point0, const XYZ& vector0, const XYZ& point1, const XYZ& vector1, double& param0, double& param1, XYZ& intersectPoint);

	};
}

