/*
 * Author:
 * 2023/06/29 - Yuqing Liang (BIMCoder Liang)
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
	class XYZ;
	class LNLIB_EXPORT Projection
	{
	public:
		static void PointToLine(const XYZ& origin, const XYZ& vector, const XYZ& Point, XYZ& ProjectPoint);
	};
}


