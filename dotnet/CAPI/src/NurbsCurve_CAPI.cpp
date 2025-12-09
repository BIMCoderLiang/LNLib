/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "NurbsCurve_CAPI.h"
#include "NurbsCurve.h"
#include "LNLibDefinitions.h"
#include "Matrix4d.h"
#include <vector>
#include <thread>

using namespace LNLib;

LN_NurbsCurve_C LNLIB_NURBSCUR_create_line(XYZ_C start, XYZ_C end) {
	LN_NurbsCurve result;
	NurbsCurve::CreateLine(ToXYZ(start), ToXYZ(end), result);
	return ToCAPI(result);
}

void LNLIB_NURBSCUR_create_cubic_hermite(const XYZ_C* through_points, int through_points_count, const XYZ_C* tangents, int tangents_count, LN_NurbsCurve_C* out_curve)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);

	std::vector<XYZ> tgts(tangents_count);
	for (int i = 0; i < tangents_count; ++i) tgts[i] = ToXYZ(tangents[i]);
	LN_NurbsCurve result;
	NurbsCurve::CreateCubicHermite(thps, tgts, result);
	*out_curve = ToCAPI(result);
}

int LNLIB_NURBSCUR_create_arc(
	XYZ_C center, XYZ_C x_axis, XYZ_C y_axis,
	double start_rad, double end_rad,
	double x_radius, double y_radius,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	bool ok = NurbsCurve::CreateArc(ToXYZ(center), ToXYZ(x_axis), ToXYZ(y_axis),
		start_rad, end_rad, x_radius, y_radius, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

bool LNLIB_NURBSCUR_create_one_conic_arc(XYZ_C start, XYZ_C start_tangent, XYZ_C end, XYZ_C end_tangent, XYZ_C point_on_conic, XYZ_C* out_project_point, double* out_project_point_weight)
{
	XYZ point;
	double weight;
	bool ok = NurbsCurve::CreateOneConicArc(ToXYZ(start), ToXYZ(start_tangent), ToXYZ(end), ToXYZ(end_tangent), ToXYZ(point_on_conic), point, weight);
	*out_project_point = FromXYZ(point);
	*out_project_point_weight = weight;
	return ok;
}

void LNLIB_NURBSCUR_split_arc(XYZ_C start, XYZ_C project_point, double project_point_weight, XYZ_C end, XYZ_C* out_insert_point_at_start_side, XYZ_C* out_split_point, XYZ_C* out_insert_point_at_end_side, double* out_insert_weight)
{
	XYZ insert_sp;
	XYZ sp;
	XYZ insert_ep;
	double weight = 0.0;
	NurbsCurve::SplitArc(ToXYZ(start), ToXYZ(project_point), project_point_weight, ToXYZ(end), insert_sp, sp, insert_ep, weight);
	*out_insert_point_at_start_side = FromXYZ(insert_sp);
	*out_split_point = FromXYZ(sp);
	*out_insert_point_at_end_side = FromXYZ(insert_ep);
	*out_insert_weight = weight;
}

int LNLIB_NURBSCUR_create_open_conic(
	XYZ_C start, XYZ_C start_tangent,
	XYZ_C end, XYZ_C end_tangent,
	XYZ_C point_on_conic,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	bool ok = NurbsCurve::CreateOpenConic(ToXYZ(start), ToXYZ(start_tangent),
		ToXYZ(end), ToXYZ(end_tangent), ToXYZ(point_on_conic), result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSCUR_global_interpolation(
	int degree,
	const XYZ_C* points,
	int point_count,
	LN_NurbsCurve_C* out_curve) {
	std::vector<XYZ> pts(point_count);
	for (int i = 0; i < point_count; ++i) pts[i] = ToXYZ(points[i]);
	LN_NurbsCurve result;
	NurbsCurve::GlobalInterpolation(degree, pts, result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_global_interpolation_with_tangents(
	int degree,
	const XYZ_C* points,
	const XYZ_C* tangents,
	double tangent_factor,
	int point_count,
	LN_NurbsCurve_C* out_curve) {
	std::vector<XYZ> pts(point_count), tans(point_count);
	for (int i = 0; i < point_count; ++i) {
		pts[i] = ToXYZ(points[i]);
		tans[i] = ToXYZ(tangents[i]);
	}
	LN_NurbsCurve result;
	NurbsCurve::GlobalInterpolation(degree, pts, tans, tangent_factor, result);
	*out_curve = ToCAPI(result);
}

int LNLIB_NURBSCUR_cubic_local_interpolation(
	const XYZ_C* points,
	int point_count,
	LN_NurbsCurve_C* out_curve) {
	std::vector<XYZ> pts(point_count);
	for (int i = 0; i < point_count; ++i) pts[i] = ToXYZ(points[i]);
	LN_NurbsCurve result;
	bool ok = NurbsCurve::CubicLocalInterpolation(pts, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

int LNLIB_NURBSCUR_least_squares_approximation(
	int degree,
	const XYZ_C* points,
	int point_count,
	int control_point_count,
	LN_NurbsCurve_C* out_curve) {
	std::vector<XYZ> pts(point_count);
	for (int i = 0; i < point_count; ++i) pts[i] = ToXYZ(points[i]);
	LN_NurbsCurve result;
	bool ok = NurbsCurve::LeastSquaresApproximation(degree, pts, control_point_count, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

int LNLIB_NURBSCUR_weighted_constrained_least_squares(
	int degree,
	const XYZ_C* points,
	const double* point_weights,
	const XYZ_C* tangents,
	const int* tangent_indices,
	const double* tangent_weights,
	int tangent_count,
	int control_point_count,
	LN_NurbsCurve_C* out_curve) {
	std::vector<XYZ> pts, tans;
	std::vector<double> pws, tws;
	std::vector<int> tidx;
	for (int i = 0; i < control_point_count; ++i) {
		pts.push_back(ToXYZ(points[i]));
		pws.push_back(point_weights[i]);
	}
	for (int i = 0; i < tangent_count; ++i) {
		tans.push_back(ToXYZ(tangents[i]));
		tidx.push_back(tangent_indices[i]);
		tws.push_back(tangent_weights[i]);
	}
	LN_NurbsCurve result;
	bool ok = NurbsCurve::WeightedAndContrainedLeastSquaresApproximation(
		degree, pts, pws, tans, tidx, tws, control_point_count, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

double LNLIB_NURBSCUR_computer_remove_knot_error_bound(LN_NurbsCurve_C curve, int removal_index)
{
	return NurbsCurve::ComputerRemoveKnotErrorBound(FromCAPI(curve), removal_index);
}

void LNLIB_NURBSCUR_remove_knots_by_given_bound(LN_NurbsCurve_C curve, const double* params, int params_count, double maxError, double* out_errors, int* out_errors_count, LN_NurbsCurve_C* out_curve)
{
	std::vector<double> ps(params_count);
	for (int i = 0; i < params_count; ++i) ps[i] = params[i];
	std::vector<double> errs;
	LN_NurbsCurve result;
	NurbsCurve::RemoveKnotsByGivenBound(FromCAPI(curve), ps, errs, maxError, result);
	for (size_t i = 0; i < errs.size(); ++i) {
		out_errors[i] = errs[i];
	}
	*out_errors_count = static_cast<int>(errs.size());
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_global_approximation_by_error_bound(
	int degree,
	const XYZ_C* points,
	int point_count,
	double max_error,
	LN_NurbsCurve_C* out_curve) {
	std::vector<XYZ> pts(point_count);
	for (int i = 0; i < point_count; ++i) pts[i] = ToXYZ(points[i]);
	LN_NurbsCurve result;
	NurbsCurve::GlobalApproximationByErrorBound(degree, pts, max_error, result);
	*out_curve = ToCAPI(result);
}

XYZ_C LNLIB_NURBSCUR_get_point_on_curve(LN_NurbsCurve_C curve, double param_t) {
	return FromXYZ(NurbsCurve::GetPointOnCurve(FromCAPI(curve), param_t));
}

XYZ_C LNLIB_NURBSCUR_get_point_on_curve_by_corner_cut(LN_NurbsCurve_C curve, double param_t) {
	return FromXYZ(NurbsCurve::GetPointOnCurveByCornerCut(FromCAPI(curve), param_t));
}

int LNLIB_NURBSCUR_compute_rational_curve_derivatives(
	LN_NurbsCurve_C curve,
	int derivative_order,
	double param_t,
	XYZ_C* out_derivatives) {
	auto derivatives = NurbsCurve::ComputeRationalCurveDerivatives(FromCAPI(curve), derivative_order, param_t);
	for (size_t i = 0; i < derivatives.size(); ++i) {
		out_derivatives[i] = FromXYZ(derivatives[i]);
	}
	return static_cast<int>(derivatives.size());
}

double LNLIB_NURBSCUR_curvature(LN_NurbsCurve_C curve, double param_t) {
	return NurbsCurve::Curvature(FromCAPI(curve), param_t);
}

double LNLIB_NURBSCUR_torsion(LN_NurbsCurve_C curve, double param_t) {
	return NurbsCurve::Torsion(FromCAPI(curve), param_t);
}

int LNLIB_NURBSCUR_fit_with_cubic(const XYZ_C* through_points, int through_points_count, int start_point_index, int end_point_index, const XYZ_C start_tangent, const XYZ_C end_tangent, double max_error, XYZW_C* middle_control_points, int* out_size)
{
	std::vector<XYZ> thps(through_points_count);
	for (int i = 0; i < through_points_count; ++i) thps[i] = ToXYZ(through_points[i]);

	std::vector<XYZW> mcps;
	bool ok = NurbsCurve::FitWithCubic(thps, start_point_index, end_point_index, ToXYZ(start_tangent), ToXYZ(end_tangent), max_error, mcps);
	*out_size = static_cast<int>(mcps.size());
	for (size_t i = 0; i < mcps.size(); ++i) {
		middle_control_points[i] = { mcps[i].WX(),mcps[i].WY(),mcps[i].WZ(),mcps[i].W() };
	}
	return ok ? 1 : 0;
}

XYZ_C LNLIB_NURBSCUR_normal(LN_NurbsCurve_C curve, LNLIB_ENUMS_CurveNormal_C normal_type, double param_t) {
	return FromXYZ(NurbsCurve::Normal(FromCAPI(curve), static_cast<CurveNormal>(normal_type), param_t));
}

int LNLIB_NURBSCUR_project_normal(LN_NurbsCurve_C curve, XYZ_C* out_normals) {
	auto normals = NurbsCurve::ProjectNormal(FromCAPI(curve));
	for (size_t i = 0; i < normals.size(); ++i) {
		out_normals[i] = FromXYZ(normals[i]);
	}
	return static_cast<int>(normals.size());
}

double LNLIB_NURBSCUR_get_param_on_curve(LN_NurbsCurve_C curve, XYZ_C given_point) {
	return NurbsCurve::GetParamOnCurve(FromCAPI(curve), ToXYZ(given_point));
}

double LNLIB_NURBSCUR_approximate_length(LN_NurbsCurve_C curve, LNLIB_ENUMS_IntegratorType_C integrator_type) {
	return NurbsCurve::ApproximateLength(FromCAPI(curve), static_cast<IntegratorType>(integrator_type));
}


double LNLIB_NURBSCUR_get_param_by_length(LN_NurbsCurve_C curve, double given_length, LNLIB_ENUMS_IntegratorType_C integrator_type) {
	return NurbsCurve::GetParamOnCurve(FromCAPI(curve), given_length, static_cast<IntegratorType>(integrator_type));
}

int LNLIB_NURBSCUR_get_params_by_equal_length(
	LN_NurbsCurve_C curve,
	double segment_length,
	LNLIB_ENUMS_IntegratorType_C integrator_type,
	double* out_params) {
	auto params = NurbsCurve::GetParamsOnCurve(FromCAPI(curve), segment_length,
		static_cast<IntegratorType>(integrator_type));
	for (size_t i = 0; i < params.size(); ++i) {
		out_params[i] = params[i];
	}
	return static_cast<int>(params.size());
}

int LNLIB_NURBSCUR_split_at(
	LN_NurbsCurve_C curve,
	double param_t,
	LN_NurbsCurve_C* out_left,
	LN_NurbsCurve_C* out_right) {
	LN_NurbsCurve left, right;
	bool ok = NurbsCurve::SplitAt(FromCAPI(curve), param_t, left, right);
	if (ok) {
		*out_left = ToCAPI(left);
		*out_right = ToCAPI(right);
	}
	return ok ? 1 : 0;
}

int LNLIB_NURBSCUR_segment(
	LN_NurbsCurve_C curve,
	double start_param,
	double end_param,
	LN_NurbsCurve_C* out_segment) {
	LN_NurbsCurve seg;
	bool ok = NurbsCurve::Segment(FromCAPI(curve), start_param, end_param, seg);
	if (ok) *out_segment = ToCAPI(seg);
	return ok ? 1 : 0;
}

int LNLIB_NURBSCUR_merge(LN_NurbsCurve_C left, LN_NurbsCurve_C right, LN_NurbsCurve_C* out_curve)
{
	LNLib::LN_NurbsCurve result;
	bool ok = NurbsCurve::Merge(FromCAPI(left), FromCAPI(right), result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSCUR_offset(LN_NurbsCurve_C curve, double offset, LNLIB_ENUMS_OffsetType_C type, LN_NurbsCurve_C* out_curve)
{
	LNLib::LN_NurbsCurve result;
	NurbsCurve::Offset(FromCAPI(curve), offset, static_cast<OffsetType>(type), result);
	*out_curve = ToCAPI(result);
}

int LNLIB_NURBSCUR_decompose_to_beziers(
	LN_NurbsCurve_C curve,
	LN_NurbsCurve_C* out_segments,
	int max_segments) {
	auto segments = NurbsCurve::DecomposeToBeziers(FromCAPI(curve));
	int n = std::min(static_cast<int>(segments.size()), max_segments);
	for (int i = 0; i < n; ++i) {
		out_segments[i] = ToCAPI(segments[i]);
	}
	return n;
}

int LNLIB_NURBSCUR_tessellate(
	LN_NurbsCurve_C curve,
	XYZ_C* out_points) {
	auto points = NurbsCurve::Tessellate(FromCAPI(curve));
	for (size_t i = 0; i < points.size(); ++i) {
		out_points[i] = FromXYZ(points[i]);
	}
	return static_cast<int>(points.size());
}

void LNLIB_NURBSCUR_create_transformed(
	LN_NurbsCurve_C curve,
	Matrix4d_C matrix,
	LN_NurbsCurve_C* out_curve) {
	const double* m = matrix.m;
	Matrix4d mat(
		m[0], m[1], m[2], m[3],
		m[4], m[5], m[6], m[7],
		m[8], m[9], m[10], m[11],
		m[12], m[13], m[14], m[15]
	);
	LN_NurbsCurve result;
	NurbsCurve::CreateTransformed(FromCAPI(curve), mat, result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_reverse(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::Reverse(FromCAPI(curve), result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_reparametrize(
	LN_NurbsCurve_C curve,
	double min_val,
	double max_val,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::Reparametrize(FromCAPI(curve), min_val, max_val, result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_reparametrize_linear_rational(
	LN_NurbsCurve_C curve,
	double alpha, double beta,
	double gamma, double delta,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::Reparametrize(FromCAPI(curve), alpha, beta, gamma, delta, result);
	*out_curve = ToCAPI(result);
}

int LNLIB_NURBSCUR_insert_knot(
	LN_NurbsCurve_C curve,
	double knot,
	int times,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	int r = NurbsCurve::InsertKnot(FromCAPI(curve), knot, times, result);
	if (r >= 0) *out_curve = ToCAPI(result);
	return r;
}

int LNLIB_NURBSCUR_remove_knot(
	LN_NurbsCurve_C curve,
	double knot,
	int times,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	bool ok = NurbsCurve::RemoveKnot(FromCAPI(curve), knot, times, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSCUR_remove_excessive_knots(
	LN_NurbsCurve_C curve,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::RemoveExcessiveKnots(FromCAPI(curve), result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_refine_knot_vector(
	LN_NurbsCurve_C curve,
	const double* insert_knots,
	int insert_count,
	LN_NurbsCurve_C* out_curve) {
	std::vector<double> knots(insert_knots, insert_knots + insert_count);
	LN_NurbsCurve result;
	NurbsCurve::RefineKnotVector(FromCAPI(curve), knots, result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_elevate_degree(
	LN_NurbsCurve_C curve,
	int times,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::ElevateDegree(FromCAPI(curve), times, result);
	*out_curve = ToCAPI(result);
}

int LNLIB_NURBSCUR_reduce_degree(
	LN_NurbsCurve_C curve,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	bool ok = NurbsCurve::ReduceDegree(FromCAPI(curve), result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

int LNLIB_NURBSCUR_is_closed(LN_NurbsCurve_C curve) {
	return NurbsCurve::IsClosed(FromCAPI(curve)) ? 1 : 0;
}

int LNLIB_NURBSCUR_is_linear(LN_NurbsCurve_C curve) {
	return NurbsCurve::IsLinear(FromCAPI(curve)) ? 1 : 0;
}

int LNLIB_NURBSCUR_is_arc(LN_NurbsCurve_C curve, LNLIB_ArcInfo_C* out_arc_info)
{
	LN_ArcInfo info;
	bool result = NurbsCurve::IsArc(FromCAPI(curve), info);
	LNLIB_ArcInfo_C r;
	r.radius = info.Radius;
	r.center = FromXYZ(info.Center);
	*out_arc_info = r;
	return result;
}

int LNLIB_NURBSCUR_is_clamp(LN_NurbsCurve_C curve) {
	return NurbsCurve::IsClamp(FromCAPI(curve)) ? 1 : 0;
}

int LNLIB_NURBSCUR_is_periodic(LN_NurbsCurve_C curve) {
	return NurbsCurve::IsPeriodic(FromCAPI(curve)) ? 1 : 0;
}

int LNLIB_NURBSCUR_can_compute_derivative(LN_NurbsCurve_C curve, double param_t) {
	return NurbsCurve::CanComputeDerivative(FromCAPI(curve), param_t) ? 1 : 0;
}

int LNLIB_NURBSCUR_control_point_reposition(
	LN_NurbsCurve_C curve,
	double param_t,
	int move_index,
	XYZ_C move_direction,
	double move_distance,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	bool ok = NurbsCurve::ControlPointReposition(FromCAPI(curve), param_t, move_index,
		ToXYZ(move_direction), move_distance, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSCUR_weight_modification(
	LN_NurbsCurve_C curve,
	double param_t,
	int move_index,
	double move_distance,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::WeightModification(FromCAPI(curve), param_t, move_index, move_distance, result);
	*out_curve = ToCAPI(result);
}

int LNLIB_NURBSCUR_neighbor_weights_modification(
	LN_NurbsCurve_C curve,
	double param_t,
	int move_index,
	double move_distance,
	double scale,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	bool ok = NurbsCurve::NeighborWeightsModification(FromCAPI(curve), param_t, move_index,
		move_distance, scale, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSCUR_warping(
	LN_NurbsCurve_C curve,
	const double* warp_shape,
	int warp_shape_count,
	double warp_distance,
	XYZ_C plane_normal,
	double start_param,
	double end_param,
	LN_NurbsCurve_C* out_curve) {
	std::vector<double> shape(warp_shape, warp_shape + warp_shape_count);
	LN_NurbsCurve result;
	NurbsCurve::Warping(FromCAPI(curve), shape, warp_distance, ToXYZ(plane_normal), start_param, end_param, result);
	*out_curve = ToCAPI(result);
}

int LNLIB_NURBSCUR_flattening(
	LN_NurbsCurve_C curve,
	XYZ_C line_start,
	XYZ_C line_end,
	double start_param,
	double end_param,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	bool ok = NurbsCurve::Flattening(FromCAPI(curve), ToXYZ(line_start), ToXYZ(line_end),
		start_param, end_param, result);
	if (ok) *out_curve = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSCUR_bending(
	LN_NurbsCurve_C curve,
	double start_param,
	double end_param,
	XYZ_C bend_center,
	double radius,
	double cross_ratio,
	LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::Bending(FromCAPI(curve), start_param, end_param, ToXYZ(bend_center), radius, cross_ratio, result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_constraint_based_modification(
	LN_NurbsCurve_C curve,
	const double* constraint_params,
	const XYZ_C* derivative_constraints,
	const int* applied_indices,
	const int* applied_degrees,
	const int* fixed_cp_indices,
	int constraint_count,
	int fixed_count,
	LN_NurbsCurve_C* out_curve) {
	std::vector<double> params(constraint_params, constraint_params + constraint_count);
	std::vector<XYZ> derivs(constraint_count);
	for (int i = 0; i < constraint_count; ++i) derivs[i] = ToXYZ(derivative_constraints[i]);
	std::vector<int> aidx(applied_indices, applied_indices + constraint_count);
	std::vector<int> adegs(applied_degrees, applied_degrees + constraint_count);
	std::vector<int> fidx(fixed_cp_indices, fixed_cp_indices + fixed_count);
	LN_NurbsCurve result;
	NurbsCurve::ConstraintBasedModification(FromCAPI(curve), params, derivs, aidx, adegs, fidx, result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_to_clamp_curve(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::ToClampCurve(FromCAPI(curve), result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_to_unclamp_curve(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
	LN_NurbsCurve result;
	NurbsCurve::ToUnclampCurve(FromCAPI(curve), result);
	*out_curve = ToCAPI(result);
}

void LNLIB_NURBSCUR_equally_tessellate(
	LN_NurbsCurve_C curve,
	int max_count,
	XYZ_C* out_points,
	double* out_knots) {
	std::vector<XYZ> pts;
	std::vector<double> knots;
	NurbsCurve::EquallyTessellate(FromCAPI(curve), pts, knots);
	int n = std::min(static_cast<int>(pts.size()), max_count);
	for (int i = 0; i < n; ++i) {
		out_points[i] = FromXYZ(pts[i]);
		out_knots[i] = knots[i];
	}
}

