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
#include "XYZ_CAPI.h"

LNLIB_ENUMS_CurveCurveIntersectionType_C LNLIB_INTERSECT_compute_rays(
	XYZ_C point_0, XYZ_C vector_0,
	XYZ_C point_1, XYZ_C vector_1,
	double* param_0,
	double* param_1,
	XYZ_C* intersect_point)
{
	double t0, t1;
	LNLib::XYZ ip;
	auto type = LNLib::Intersection::ComputeRays(
		ToXYZ(point_0), ToXYZ(vector_0),
		ToXYZ(point_1), ToXYZ(vector_1),
		t0, t1, ip);

	if (param_0) *param_0 = t0;
	if (param_1) *param_1 = t1;
	if (intersect_point) *intersect_point = FromXYZ(ip);

	return static_cast<LNLIB_ENUMS_CurveCurveIntersectionType_C>(static_cast<int>(type));
}

LNLIB_ENUMS_LinePlaneIntersectionType_C LNLIB_INTERSECT_compute_line_and_plane(
	XYZ_C normal,
	XYZ_C point_on_plane,
	XYZ_C point_on_line,
	XYZ_C line_direction,
	XYZ_C* intersect_point)
{
	LNLib::XYZ ip;
	auto type = LNLib::Intersection::ComputeLineAndPlane(
		ToXYZ(normal), ToXYZ(point_on_plane),
		ToXYZ(point_on_line), ToXYZ(line_direction), ip);

	if (intersect_point) *intersect_point = FromXYZ(ip);
	return static_cast<LNLIB_ENUMS_LinePlaneIntersectionType_C>(static_cast<int>(type));
}