/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Intersection_CAPI.h"
#include "Intersection.h"
#include "LNEnums.h"
#include "XYZ.h"

extern "C" {

	static LNLib::XYZ ToXYZ(XYZ_C c) { return LNLib::XYZ(c.x, c.y, c.z); }
	static XYZ_C FromXYZ(const LNLib::XYZ& v) { return { v.GetX(), v.GetY(), v.GetZ() }; }

	LNLIB_EXPORT CurveCurveIntersectionType_C intersection_compute_rays(
		XYZ_C point0, XYZ_C vector0,
		XYZ_C point1, XYZ_C vector1,
		double* out_param0,
		double* out_param1,
		XYZ_C* out_intersect_point)
	{
		double t0, t1;
		LNLib::XYZ ip;
		auto type = LNLib::Intersection::ComputeRays(
			ToXYZ(point0), ToXYZ(vector0),
			ToXYZ(point1), ToXYZ(vector1),
			t0, t1, ip);

		if (out_param0) *out_param0 = t0;
		if (out_param1) *out_param1 = t1;
		if (out_intersect_point) *out_intersect_point = FromXYZ(ip);

		return static_cast<CurveCurveIntersectionType_C>(static_cast<int>(type));
	}

	LNLIB_EXPORT LinePlaneIntersectionType_C intersection_compute_line_and_plane(
		XYZ_C plane_normal,
		XYZ_C point_on_plane,
		XYZ_C point_on_line,
		XYZ_C line_direction,
		XYZ_C* out_intersect_point)
	{
		LNLib::XYZ ip;
		auto type = LNLib::Intersection::ComputeLineAndPlane(
			ToXYZ(plane_normal), ToXYZ(point_on_plane),
			ToXYZ(point_on_line), ToXYZ(line_direction), ip);

		if (out_intersect_point) *out_intersect_point = FromXYZ(ip);
		return static_cast<LinePlaneIntersectionType_C>(static_cast<int>(type));
	}

}