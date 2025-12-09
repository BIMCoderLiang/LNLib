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

int LNLIB_KV_get_continuity(
	int degree,
	const double* knot_vector,
	int knot_vector_count,
	double knot)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	return LNLib::KnotVectorUtils::GetContinuity(degree, kv, knot);
}

void LNLIB_KV_rescale(
	const double* knot_vector,
	int knot_vector_count,
	double min,
	double max,
	double* result)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	auto res = LNLib::KnotVectorUtils::Rescale(kv, min, max);
	std::memcpy(result, res.data(), knot_vector_count * sizeof(double));
}

int LNLIB_KV_is_uniform(
	const double* knot_vector,
	int knot_vector_count)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	return LNLib::KnotVectorUtils::IsUniform(kv) ? 1 : 0;
}

void LNLIB_KV_get_knot_multiplicity_map(
	const double* knot_vector,
	int knot_vector_count,
	int* out_size,
	double* out_keys,
	int* out_values)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	auto map = LNLib::KnotVectorUtils::GetKnotMultiplicityMap(kv);
	*out_size = static_cast<int>(map.size());
	int i = 0;
	for (const auto& pair : map) {
		out_keys[i] = pair.first;
		out_values[i] = pair.second;
		++i;
	}
}

void LNLIB_KV_get_internal_knot_multiplicity_map(
	const double* knot_vector,
	int knot_vector_count,
	int* out_size,
	double* out_keys,
	int* out_values)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	auto map = LNLib::KnotVectorUtils::GetInternalKnotMultiplicityMap(kv);
	*out_size = static_cast<int>(map.size());
	int i = 0;
	for (const auto& pair : map) {
		out_keys[i] = pair.first;
		out_values[i] = pair.second;
		++i;
	}
}