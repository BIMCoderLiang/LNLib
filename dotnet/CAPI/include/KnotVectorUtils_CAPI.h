/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

LNLIB_EXPORT int LNLIB_KV_get_continuity(
    int degree,
    const double* knot_vector,
    int knot_vector_count,
    double knot);

LNLIB_EXPORT void LNLIB_KV_rescale(
    const double* knot_vector,
    int knot_vector_count,
    double min,
    double max,
    double* result);

LNLIB_EXPORT int LNLIB_KV_is_uniform(
    const double* knot_vector,
    int knot_vector_count);

LNLIB_EXPORT void LNLIB_KV_get_knot_multiplicity_map(
    const double* knot_vector,
    int knot_vector_count,
    int* out_size,
    double* out_keys,
    int* out_values);

LNLIB_EXPORT void LNLIB_KV_get_internal_knot_multiplicity_map(
    const double* knot_vector,
    int knot_vector_count,
    int* out_size,
    double* out_keys,
    int* out_values);

#ifdef __cplusplus
}
#endif