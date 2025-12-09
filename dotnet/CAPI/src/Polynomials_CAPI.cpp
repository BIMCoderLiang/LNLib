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
#include "Constants.h"
#include <vector>
#include <cstring>


double LNLIB_POLY_bernstein(int index, int degree, double param_t)
{
	return LNLib::Polynomials::Bernstein(index, degree, param_t);
}

void LNLIB_POLY_all_bernstein(int degree, double param_t, double* out_array)
{
	auto result = LNLib::Polynomials::AllBernstein(degree, param_t);
	std::memcpy(out_array, result.data(), (degree + 1) * sizeof(double));
}

void LNLIB_POLY_horner_uv(
	int degree_u,
	int degree_v,
	const double* coefficients,
	const UV_C* uv,
	double* out_horner)
{
	std::vector<std::vector<double>> coeffs(degree_u + 1, std::vector<double>(degree_v + 1));
	for (int i = 0; i <= degree_u; ++i) {
		for (int j = 0; j <= degree_v; ++j) {
			coeffs[i][j] = coefficients[i * (degree_v + 1) + j];
		}
	}

	LNLib::UV cpp_uv{ uv->u, uv->v };
	double val = LNLib::Polynomials::Horner(degree_u, degree_v, coeffs, cpp_uv);
	*out_horner = val;
}

double LNLIB_POLY_horner(int degree, const double* coefficients, int coeff_count, double param_t)
{
	std::vector<double> coeffs(coefficients, coefficients + coeff_count);
	return LNLib::Polynomials::Horner(degree, coeffs, param_t);
}

int LNLIB_POLY_get_knot_span_index(int degree, const double* knot_vector, int knot_vector_count, double param_t)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	return LNLib::Polynomials::GetKnotSpanIndex(degree, kv, param_t);
}

int LNLIB_POLY_get_knot_multiplicity(const double* knot_vector, int knot_vector_count, double param_t)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	return LNLib::Polynomials::GetKnotMultiplicity(kv, param_t);
}

void LNLIB_POLY_basis_functions(
	int span_index,
	int degree,
	const double* knot_vector,
	int knot_vector_count,
	double param_t,
	double* out_basis_functions
)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	LNLib::Polynomials::BasisFunctions(span_index, degree, kv, param_t, out_basis_functions);
}

void LNLIB_POLY_basis_functions_derivatives(
	int span_index,
	int degree,
	int derivative,
	const double* knot_vector,
	int knot_vector_count,
	double param_t,
	double** out_derivatives,
	int* out_rows,
	int* out_cols
)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	auto result = LNLib::Polynomials::BasisFunctionsDerivatives(span_index, degree, derivative, kv, param_t);

	if (result.empty()) {
		*out_derivatives = nullptr;
		*out_rows = 0;
		*out_cols = 0;
		return;
	}

	int rows = static_cast<int>(result.size());
	int cols = rows > 0 ? static_cast<int>(result[0].size()) : 0;

	double* flat = new(std::nothrow) double[rows * cols];
	if (!flat) {
		return;
	}

	for (int i = 0; i < rows; ++i) {
		if (static_cast<int>(result[i].size()) != cols) {
			delete[] flat;
			return;
		}
		for (int j = 0; j < cols; ++j) {
			flat[i * cols + j] = result[i][j];
		}
	}

	*out_derivatives = flat;
	*out_rows = rows;
	*out_cols = cols;
}

void LNLIB_POLY_basis_functions_first_derivative(
	int span_index,
	int degree,
	const double* knot_vector,
	int knot_vector_count,
	double param_t,
	double* out_derivatives
)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	double der[2][LNLib::Constants::NURBSMaxDegree + 1];
	LNLib::Polynomials::BasisFunctionsFirstOrderDerivative(span_index, degree, kv, param_t, der);

	int cols = degree + 1;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < cols; ++j) {
			out_derivatives[i * (LNLib::Constants::NURBSMaxDegree + 1) + j] = der[i][j];
		}
		for (int j = cols; j <= LNLib::Constants::NURBSMaxDegree; ++j) {
			out_derivatives[i * (LNLib::Constants::NURBSMaxDegree + 1) + j] = 0.0;
		}
	}
}

double LNLIB_POLY_one_basis_function(int span_index, int degree, const double* knot_vector, int knot_vector_count, double param_t)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	double val = LNLib::Polynomials::OneBasisFunction(span_index, degree, kv, param_t);
	return val;
}

void LNLIB_POLY_one_basis_function_derivative(
	int span_index,
	int degree,
	int derivative,
	const double* knot_vector,
	int knot_vector_count,
	double param_t,
	double** out_derivatives,
	int* out_size)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);
	std::vector<double> result = LNLib::Polynomials::OneBasisFunctionDerivative(span_index, degree, derivative, kv, param_t);

	int size = static_cast<int>(result.size());
	if (size != derivative + 1) {
		return;
	}
	double* data = new(std::nothrow) double[size];
	if (!data) {
		return;
	}

	for (int i = 0; i < size; ++i) {
		data[i] = result[i];
	}

	*out_derivatives = data;
	*out_size = size;
}

void LNLIB_POLY_all_basis_functions(
	int span_index,
	int degree,
	const double* knot_vector,
	int knot_vector_count,
	double knot,
	double** out_basis_functions,
	int* out_rows,
	int* out_cols
)
{
	std::vector<double> kv(knot_vector, knot_vector + knot_vector_count);

	auto result = LNLib::Polynomials::AllBasisFunctions(span_index, degree, kv, knot);

	if (result.empty()) {
		*out_basis_functions = nullptr;
		*out_rows = 0;
		*out_cols = 0;
		return;
	}

	int rows = static_cast<int>(result.size());
	int cols = rows > 0 ? static_cast<int>(result[0].size()) : 0;

	for (int i = 0; i < rows; ++i) {
		if (static_cast<int>(result[i].size()) != cols) {
			return;
		}
	}

	double* flat = new(std::nothrow) double[rows * cols];
	if (!flat) {
		return;
	}

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			flat[i * cols + j] = result[i][j];
		}
	}

	*out_basis_functions = flat;
	*out_rows = rows;
	*out_cols = cols;
}

void LNLIB_POLY_bezier_to_power_matrix(
	int degree,
	double** out_matrix,
	int* out_rows,
	int* out_cols)
{
	auto matrix = LNLib::Polynomials::BezierToPowerMatrix(degree);

	int n = static_cast<int>(matrix.size());
	if (n == 0) {
		*out_matrix = nullptr;
		*out_rows = 0;
		*out_cols = 0;
		return;
	}
	if (n != degree + 1) {
		return;
	}

	int cols = n;
	for (int i = 0; i < n; ++i) {
		if (static_cast<int>(matrix[i].size()) != cols) {
			return;
		}
	}

	double* flat = new(std::nothrow) double[n * cols];
	if (!flat) {
		return;
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < cols; ++j) {
			flat[i * cols + j] = matrix[i][j];
		}
	}

	*out_matrix = flat;
	*out_rows = n;
	*out_cols = cols;
}

void LNLIB_POLY_power_to_bezier_matrix(
	int degree,
	const double* matrix,
	double** out_matrix,
	int* out_rows,
	int* out_cols)
{
	const int n = degree + 1;
	const int total_size = n * n;

	std::vector<std::vector<double>> input(n, std::vector<double>(n));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			input[i][j] = matrix[i * n + j];
		}
	}

	auto result = LNLib::Polynomials::PowerToBezierMatrix(degree, input);
	if (static_cast<int>(result.size()) != n) {
		return;
	}
	for (int i = 0; i < n; ++i) {
		if (static_cast<int>(result[i].size()) != n) {
			return;
		}
	}

	double* flat = new(std::nothrow) double[total_size];
	if (!flat) {
		return;
	}

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			flat[i * n + j] = out_matrix[i][j];
		}
	}

	*out_matrix = flat;
	*out_rows = n;
	*out_cols = n;
}

