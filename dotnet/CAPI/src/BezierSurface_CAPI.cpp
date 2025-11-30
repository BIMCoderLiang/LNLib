/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "BezierSurface_CAPI.h"
#include "BezierSurface.h"
#include "LNObject.h"
#include "XYZW_CAPI.h"
#include <vector>

extern "C" {

	static LNLib::XYZ ToXYZ(XYZ_C c) { return LNLib::XYZ(c.x, c.y, c.z); }
	static XYZ_C FromXYZ(const LNLib::XYZ& v) { return { v.GetX(), v.GetY(), v.GetZ() }; }

	static LNLib::XYZW ToXYZW(XYZW_C c) { return LNLib::XYZW(c.wx, c.wy, c.wz, c.w); }
	static XYZW_C FromXYZW(const LNLib::XYZW& v) { return { v.GetWX(), v.GetWY(), v.GetWZ(), v.GetW() }; }

	LNLIB_EXPORT XYZ_C bezier_surface_get_point_by_de_casteljau(
		int degree_u, int degree_v,
		const XYZ_C* control_points,
		int rows, int cols,
		UV_C uv)
	{
		std::vector<std::vector<LNLib::XYZ>> cps(rows, std::vector<LNLib::XYZ>(cols));
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				cps[i][j] = ToXYZ(control_points[i * cols + j]);

		LNLib::LN_BezierSurface<LNLib::XYZ> surf{ degree_u, degree_v, cps };
		LNLib::XYZ pt = LNLib::BezierSurface::GetPointOnSurfaceByDeCasteljau(surf, LNLib::UV(uv.u, uv.v));
		return FromXYZ(pt);
	}

}