/*
 * Author:
 * 2023/11/15 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "XYZW.h"
#include <vector>

namespace LNLib
{
	struct LNLIB_EXPORT LN_Curve
	{
		int Degree;
		std::vector<double> KnotVector;
		std::vector<XYZW> ControlPoints;
	};

	struct LNLIB_EXPORT LN_Surface
	{
		int DegreeU;
		int DegreeV;
		std::vector<double> KnotVectorU;
		std::vector<double> KnotVectorV;
		std::vector<std::vector<XYZW>> ControlPoints;
	};
}



