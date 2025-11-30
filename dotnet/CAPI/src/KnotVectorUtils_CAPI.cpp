/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "KnotVectorUtils_CAPI.h"
#include "KnotVectorUtils.h"
#include <vector>
#include <map>

extern "C" {

	LNLIB_EXPORT int knot_vector_utils_get_continuity(
		int degree,
		const double* knot_vector,
		int knot_count,
		double knot)
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		return LNLib::KnotVectorUtils::GetContinuity(degree, kv, knot);
	}

	LNLIB_EXPORT void knot_vector_utils_rescale(
		const double* knot_vector,
		int knot_count,
		double min_val,
		double max_val,
		double* out_rescaled)
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		auto res = LNLib::KnotVectorUtils::Rescale(kv, min_val, max_val);
		std::memcpy(out_rescaled, res.data(), knot_count * sizeof(double));
	}

	LNLIB_EXPORT int knot_vector_utils_is_uniform(const double* knot_vector, int knot_count)
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		return LNLib::KnotVectorUtils::IsUniform(kv) ? 1 : 0;
	}

	LNLIB_EXPORT int knot_vector_utils_get_knot_multiplicity_map_size(const double* knot_vector, int knot_count)
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		auto map = LNLib::KnotVectorUtils::GetKnotMultiplicityMap(kv);
		return static_cast<int>(map.size());
	}

	LNLIB_EXPORT void knot_vector_utils_get_knot_multiplicity_map(
		const double* knot_vector, int knot_count,
		double* out_knots,   // size = N
		int* out_multiplicities) // size = N
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		auto map = LNLib::KnotVectorUtils::GetKnotMultiplicityMap(kv);
		int i = 0;
		for (const auto& pair : map) {
			out_knots[i] = pair.first;
			out_multiplicities[i] = pair.second;
			++i;
		}
	}

}