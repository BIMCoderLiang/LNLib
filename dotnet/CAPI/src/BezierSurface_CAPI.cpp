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
#include "XYZ_CAPI.h"
#include "XYZW_CAPI.h"
#include <vector>

XYZ_C LNLIB_BEZIERSURF_get_point_on_surface_by_deCasteljau(
	int degree_u,
	int degree_v,
	const XYZ_C* control_points,
	int cp_rows,
	int cp_cols,
	UV_C uv)
{
	std::vector<std::vector<LNLib::XYZ>> cps(cp_rows, std::vector<LNLib::XYZ>(cp_cols));
	for (int i = 0; i < cp_rows; ++i)
		for (int j = 0; j < cp_cols; ++j)
			cps[i][j] = ToXYZ(control_points[i * cp_cols + j]);

	LNLib::LN_BezierSurface<LNLib::XYZ> surf{ degree_u, degree_v, cps };
	LNLib::XYZ pt = LNLib::BezierSurface::GetPointOnSurfaceByDeCasteljau(surf, LNLib::UV(uv.u, uv.v));
	return FromXYZ(pt);
}

XYZW_C LNLIB_BEZIERSURF_get_rational_point_on_surface_by_deCasteljau(
	int degree_u,
	int degree_v,
	const XYZW_C* control_points,
	int cp_rows,
	int cp_cols,
	UV_C uv)
{
	std::vector<std::vector<LNLib::XYZW>> cps(cp_rows, std::vector<LNLib::XYZW>(cp_cols));
	for (int i = 0; i < cp_rows; ++i)
		for (int j = 0; j < cp_cols; ++j)
			cps[i][j] = ToXYZW(control_points[i * cp_cols + j]);

	LNLib::LN_BezierSurface<LNLib::XYZW> surf{ degree_u, degree_v, cps };
	LNLib::XYZW pt = LNLib::BezierSurface::GetPointOnSurfaceByDeCasteljau(surf, LNLib::UV(uv.u, uv.v));
	return FromXYZW(pt);
}