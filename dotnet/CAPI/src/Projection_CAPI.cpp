/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Projection_CAPI.h"
#include "Projection.h"
#include "XYZ.h"

XYZ_C LNLIB_PROJ_point_to_ray(
	XYZ_C origin,
	XYZ_C direction,
	XYZ_C point
)
{
	LNLib::XYZ res = LNLib::Projection::PointToRay(ToXYZ(origin), ToXYZ(direction), ToXYZ(point));
	return FromXYZ(res);
}

int LNLIB_PROJ_point_to_line(
	XYZ_C start,
	XYZ_C end,
	XYZ_C point,
	XYZ_C* project_point
)
{
	LNLib::XYZ proj;
	bool success = LNLib::Projection::PointToLine(ToXYZ(start), ToXYZ(end), ToXYZ(point), proj);
	if (success && project_point) {
		*project_point = FromXYZ(proj);
	}
	return success ? 1 : 0;
}

XYZ_C LNLIB_PROJ_stereographic(
	XYZ_C point_on_sphere,
	double radius
)
{
	LNLib::XYZ res = LNLib::Projection::Stereographic(ToXYZ(point_on_sphere), radius);
	return FromXYZ(res);
}

