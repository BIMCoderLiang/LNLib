/*
 * Author:
 * 2023/06/08 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 微信公众号：BIMCoder梁老师
 *
 * Use of this source code is governed by a GPL-3.0 license that can be found in
 * the LICENSE file.
 */

#pragma once

#include "LNLibDefinitions.h"
#include <vector>

namespace LNLib {

	class UV;
	class XYZ;
	class XYZW;
	class LNLIB_EXPORT BezierSurface
	{
	public:

		/// <summary>
		/// The NURBS Book 2nd Edition Page39
		/// Algorithm A1.7
		/// Compute a point on a Bezier surface by the deCasteljau.
		/// 
		/// Controlpoints with (n+1) rows * (m+1) columns
		///  
		///  [0][0]  [0][1] ... ...  [0][m]     ------- v direction
		///  [1][0]  [1][1] ... ...  [1][m]    |
		///    .                               |
		///    .                               u direction
		///    .							   
		///  [n][0]  [n][1] ... ...  [n][m]      
		/// 
		/// </summary>
		static XYZ GetPointOnSurfaceByDeCasteljau(int degreeU, int degreeV, const std::vector<std::vector<XYZ>>& controlPoints, UV uv);

		/// <summary>
		/// The NURBS Book 2nd Edition Page40
		/// Compute a point on a Rational Bezier surface by the deCasteljau.
		/// </summary>
		static XYZW GetPointOnRationalSurfaceByDeCasteljau(int degreeU, int degreeV, const std::vector<std::vector<XYZW>>& controlPoints, UV uv);
	};
}



