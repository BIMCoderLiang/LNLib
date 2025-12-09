/*
 * Author:
 * 2025/11/29 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "XYZ_CAPI.h"
#include "UV_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

LNLIB_EXPORT double LNLIB_POLY_horner(int degree, const double* coefficients, int coeff_count, double param_t);
LNLIB_EXPORT double LNLIB_POLY_bernstein(int index, int degree, double param_t);
LNLIB_EXPORT void   LNLIB_POLY_all_bernstein(int degree, double param_t, double* out_array);
LNLIB_EXPORT void   LNLIB_POLY_horner_uv(
    int degree_u,
    int degree_v,
    const double* coefficients,
    const UV_C* uv,
    double* out_horner);
LNLIB_EXPORT int    LNLIB_POLY_get_knot_multiplicity(const double* knot_vector, int knot_vector_count, double param_t);
LNLIB_EXPORT int    LNLIB_POLY_get_knot_span_index(int degree, const double* knot_vector, int knot_vector_count, double param_t);
LNLIB_EXPORT void   LNLIB_POLY_basis_functions(
    int span_index, 
    int degree,
    const double* knot_vector, 
    int knot_vector_count,
    double param_t,
    double* out_basis_functions
);
LNLIB_EXPORT void   LNLIB_POLY_basis_functions_derivatives(
    int span_index,
    int degree,
    int derivative,
    const double* knot_vector,
    int knot_vector_count,
    double param_t,
    double** out_derivatives,
    int* out_rows,
    int* out_cols
);
LNLIB_EXPORT void   LNLIB_POLY_basis_functions_first_derivative(
    int span_index,
    int degree,
    const double* knot_vector,
    int knot_vector_count,
    double param_t,
    double* out_derivatives
);
LNLIB_EXPORT double LNLIB_POLY_one_basis_function(int span_index, int degree, const double* knot_vector, int knot_vector_count, double param_t);
LNLIB_EXPORT void   LNLIB_POLY_one_basis_function_derivative(
    int span_index,
    int degree,
    int derivative,
    const double* knot_vector,
    int knot_vector_count,
    double param_t,
    double** out_derivatives,
    int* out_size);
LNLIB_EXPORT void   LNLIB_POLY_all_basis_functions(
    int span_index,
    int degree,
    const double* knot_vector,
    int knot_vector_count,
    double knot,
    double** out_basis_functions,
    int* out_rows,
    int* out_cols
);
LNLIB_EXPORT void   LNLIB_POLY_bezier_to_power_matrix(
    int degree,
    double** out_matrix,
    int* out_rows,
    int* out_cols);
LNLIB_EXPORT void   LNLIB_POLY_power_to_bezier_matrix(
    int degree,
    const double* matrix,
    double** out_matrix,
    int* out_rows,
    int* out_cols);

#ifdef __cplusplus
}
#endif