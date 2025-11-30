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

extern "C" {

	static LNLib::XYZ ToXYZ(XYZ_C c) {
		return LNLib::XYZ(c.x, c.y, c.z);
	}

	static XYZ_C FromXYZ(const LNLib::XYZ& v) {
		return { v.GetX(), v.GetY(), v.GetZ() };
	}

	LNLIB_EXPORT XYZ_C projection_point_to_ray(XYZ_C origin, XYZ_C direction, XYZ_C point)
	{
		LNLib::XYZ res = LNLib::Projection::PointToRay(ToXYZ(origin), ToXYZ(direction), ToXYZ(point));
		return FromXYZ(res);
	}

	LNLIB_EXPORT int projection_point_to_line(XYZ_C start, XYZ_C end, XYZ_C point, XYZ_C* out_project_point)
	{
		LNLib::XYZ proj;
		bool success = LNLib::Projection::PointToLine(ToXYZ(start), ToXYZ(end), ToXYZ(point), proj);
		if (success && out_project_point) {
			*out_project_point = FromXYZ(proj);
		}
		return success ? 1 : 0;
	}

	LNLIB_EXPORT XYZ_C projection_stereographic(XYZ_C point_on_sphere, double radius)
	{
		LNLib::XYZ res = LNLib::Projection::Stereographic(ToXYZ(point_on_sphere), radius);
		return FromXYZ(res);
	}

}