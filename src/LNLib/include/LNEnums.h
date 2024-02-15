/*
 * Author:
 * 2024/01/13 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once

namespace LNLib
{
	enum class CurveCurveIntersectionType : int
	{
		Intersecting = 0,
		Parallel = 1,
		Coincident = 2,
		Skew = 3,
	};

	enum class LinePlaneIntersectionType : int
	{
		Intersecting = 0,
		Parallel = 1,
		On = 2,
	};

	enum class CurveNormal :int
	{
		Normal = 0,
		Binormal = 1,
	};

	enum class SurfaceDirection : int
	{
		All = 0,
		UDirection = 1,
		VDirection = 2,
	};

	enum class SurfaceCurvature : int
	{
		Maximum = 0,
		Minimum = 1,
		Gauss = 2,
		Mean = 3,
		Abs = 4,
		Rms = 5
	};

	enum class IntegratorType :int
	{
		Simpson = 0,
		GaussLegendre = 1,
		Chebyshev = 2,
	};

}



