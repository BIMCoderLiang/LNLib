/*
 * Author:
 * 2025/12/07 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "Interpolation_CAPI.h"
#include "Interpolation.h"
#include "XYZ_CAPI.h"

using namespace LNLib;

double LNLIB_INTERPOLATION_get_total_chord_length(XYZ_C* through_points, int through_points_count)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);
	return Interpolation::GetTotalChordLength(thps);
}

void LNLIB_INTERPOLATION_get_chord_parameterization(XYZ_C* through_points, int through_points_count, double* out_params, int* out_size)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);

	std::vector<double> ps = Interpolation::GetChordParameterization(thps);
	for (size_t i = 0; i < ps.size(); ++i) {
		out_params[i] = ps[i];
	}
	*out_size = static_cast<int>(ps.size());
}

double LNLIB_INTERPOLATION_get_centripetal_length(XYZ_C* through_points, int through_points_count)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);
	return Interpolation::GetCentripetalLength(thps);
}

void LNLIB_INTERPOLATION_get_centripetal_parameterization(XYZ_C* through_points, int through_points_count, double* out_params, int* out_size)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);

	std::vector<double> ps = Interpolation::GetCentripetalParameterization(thps);
	for (size_t i = 0; i < ps.size(); ++i) {
		out_params[i] = ps[i];
	}
	*out_size = static_cast<int>(ps.size());
}

void LNLIB_INTERPOLATION_average_knot_vector(int degree, double* params, int params_count, double* out_knot_vector, int* out_knot_vector_size)
{
	std::vector<double> input(params_count);
	for (int i = 0; i < params_count; i++)
	{
		input[i] = params[i];
	}

	std::vector<double> kv = Interpolation::AverageKnotVector(degree, input);
	for (size_t i = 0; i < kv.size(); ++i) {
		out_knot_vector[i] = kv[i];
	}
	*out_knot_vector_size = static_cast<int>(kv.size());
}

int LNLIB_INTERPOLATION_get_surface_mesh_parameterization(XYZ_C* through_points, int rows, int cols, double* params_u, int* out_params_u_size, double* params_v, int* out_params_v_size)
{
	std::vector<std::vector<XYZ>> pts(rows, std::vector<XYZ>(cols));
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			pts[i][j] = ToXYZ(through_points[i * cols + j]);
		}
	}
	std::vector<double> paramsU;
	std::vector<double> paramsV;
	bool ok = Interpolation::GetSurfaceMeshParameterization(pts, paramsU, paramsV);
	for (size_t i = 0; i < paramsU.size(); ++i) {
		params_u[i] = paramsU[i];
	}
	*out_params_u_size = static_cast<int>(paramsU.size());
	for (size_t i = 0; i < paramsV.size(); ++i) {
		params_v[i] = paramsV[i];
	}
	*out_params_v_size = static_cast<int>(paramsV.size());
	return ok ? 1 : 0;
}

void LNLIB_INTERPOLATION_compute_tangent_by_three_points(XYZ_C* through_points, int through_points_count, XYZ_C* out_tangents, int* out_size)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);

	std::vector<XYZ> tangents = Interpolation::ComputeTangent(thps);
	for (size_t i = 0; i < tangents.size(); ++i) {
		out_tangents[i] = FromXYZ(tangents[i]);
	}
	*out_size = static_cast<int>(tangents.size());
}

int LNLIB_INTERPOLATION_compute_tangent_by_five_points(XYZ_C* through_points, int through_points_count, XYZ_C* out_tangents, int* out_size)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);

	std::vector<XYZ> tangents;
	bool ok = Interpolation::ComputeTangent(thps, tangents);
	for (size_t i = 0; i < tangents.size(); ++i) {
		out_tangents[i] = FromXYZ(tangents[i]);
	}
	*out_size = static_cast<int>(tangents.size());
	return ok ? 1 : 0;
}

int LNLIB_INTERPOLATION_computer_weight_for_rational_quadratic_interpolation(XYZ_C start_point, XYZ_C middle_control_point, XYZ_C end_point, double* out_weight)
{
	double weight = 0.0;
	bool ok = Interpolation::ComputerWeightForRationalQuadraticInterpolation(ToXYZ(start_point), ToXYZ(middle_control_point), ToXYZ(end_point), weight);
	*out_weight = weight;
	return ok ? 1 : 0;
}

void LNLIB_INTERPOLATION_compute_knot_vector(int degree, int control_points_count, double* params, int params_count, double* out_knot_vector, int* out_knot_vector_size)
{
	std::vector<double> input(params_count);
	for (int i = 0; i < params_count; i++)
	{
		input[i] = params[i];
	}

	std::vector<double> kv = Interpolation::ComputeKnotVector(degree, control_points_count, input);
	for (size_t i = 0; i < kv.size(); ++i) {
		out_knot_vector[i] = kv[i];
	}
	*out_knot_vector_size = static_cast<int>(kv.size());
}
