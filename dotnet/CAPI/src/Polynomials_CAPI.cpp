/*
 * Author:
 * 2025/11/29 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Polynomials_CAPI.h"
#include "Polynomials.h"
#include "LNObject.h"
#include <vector>
#include <cstring>

extern "C" {

	LNLIB_EXPORT double polynomials_bernstein(int index, int degree, double paramT)
	{
		return LNLib::Polynomials::Bernstein(index, degree, paramT);
	}

	LNLIB_EXPORT void polynomials_all_bernstein(int degree, double paramT, double* out_array)
	{
		auto result = LNLib::Polynomials::AllBernstein(degree, paramT);
		std::memcpy(out_array, result.data(), (degree + 1) * sizeof(double));
	}

	LNLIB_EXPORT double polynomials_horner_curve(int degree, const double* coefficients, int coeff_count, double paramT)
	{
		std::vector<double> coeffs(coefficients, coefficients + coeff_count);
		return LNLib::Polynomials::Horner(degree, coeffs, paramT);
	}

	LNLIB_EXPORT int polynomials_get_knot_span_index(int degree, const double* knot_vector, int knot_count, double paramT)
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		return LNLib::Polynomials::GetKnotSpanIndex(degree, kv, paramT);
	}

	LNLIB_EXPORT int polynomials_get_knot_multiplicity(const double* knot_vector, int knot_count, double paramT)
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		return LNLib::Polynomials::GetKnotMultiplicity(kv, paramT);
	}

	LNLIB_EXPORT void polynomials_basis_functions(
		int span_index, int degree,
		const double* knot_vector, int knot_count,
		double paramT,
		double* basis_functions)
	{
		std::vector<double> kv(knot_vector, knot_vector + knot_count);
		LNLib::Polynomials::BasisFunctions(span_index, degree, kv, paramT, basis_functions);
	}

	LNLIB_EXPORT void polynomials_bezier_to_power_matrix(int degree, double* out_matrix)
	{
		auto matrix = LNLib::Polynomials::BezierToPowerMatrix(degree);
		int n = degree + 1;
		for (int i = 0; i < n; ++i) {
			std::memcpy(&out_matrix[i * n], matrix[i].data(), n * sizeof(double));
		}
	}

}