/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "ValidationUtils_CAPI.h"
#include "ValidationUtils.h"
#include "XYZW.h"
#include <vector>

extern "C" {

	LNLIB_EXPORT int validation_utils_is_valid_knot_vector(const double* knot_vector, int count)
	{
		std::vector<double> kv(knot_vector, knot_vector + count);
		return LNLib::ValidationUtils::IsValidKnotVector(kv) ? 1 : 0;
	}

	LNLIB_EXPORT int validation_utils_is_valid_bezier(int degree, int control_points_count)
	{
		return LNLib::ValidationUtils::IsValidBezier(degree, control_points_count) ? 1 : 0;
	}

	LNLIB_EXPORT int validation_utils_is_valid_bspline(int degree, int knot_count, int cp_count)
	{
		return LNLib::ValidationUtils::IsValidBspline(degree, knot_count, cp_count) ? 1 : 0;
	}

	LNLIB_EXPORT int validation_utils_is_valid_nurbs(int degree, int knot_count, int weighted_cp_count)
	{
		return LNLib::ValidationUtils::IsValidNurbs(degree, knot_count, weighted_cp_count) ? 1 : 0;
	}

	LNLIB_EXPORT double validation_utils_compute_curve_modify_tolerance(const XYZW_C* control_points, int count)
	{
		std::vector<LNLib::XYZW> cps;
		cps.reserve(count);
		for (int i = 0; i < count; ++i) {
			cps.emplace_back(control_points[i].wx, control_points[i].wy, control_points[i].wz, control_points[i].w);
		}
		return LNLib::ValidationUtils::ComputeCurveModifyTolerance(cps);
	}

}
