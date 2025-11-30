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

static XYZ ToCPP(XYZ_C v) { return XYZ(v.x, v.y, v.z); }
static XYZ_C ToCAPI(const XYZ& v) { return { v.X(), v.Y(), v.Z() }; }

extern "C" {

	LN_NurbsCurve_C nurbs_curve_create_line(XYZ_C start, XYZ_C end) {
		LN_NurbsCurve result;
		NurbsCurve::CreateLine(ToCPP(start), ToCPP(end), result);
		return ToCAPI(result);
	}

	int nurbs_curve_create_arc(XYZ_C center, XYZ_C x_axis, XYZ_C y_axis,
		double start_rad, double end_rad, double x_radius, double y_radius,
		LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		bool ok = NurbsCurve::CreateArc(ToCPP(center), ToCPP(x_axis), ToCPP(y_axis),
			start_rad, end_rad, x_radius, y_radius, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	int nurbs_curve_create_open_conic(XYZ_C start, XYZ_C start_tangent,
		XYZ_C end, XYZ_C end_tangent, XYZ_C point_on_conic,
		LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		bool ok = NurbsCurve::CreateOpenConic(ToCPP(start), ToCPP(start_tangent),
			ToCPP(end), ToCPP(end_tangent), ToCPP(point_on_conic), result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	void nurbs_curve_global_interpolation(int degree, const XYZ_C* points, int point_count, LN_NurbsCurve_C* out_curve) {
		std::vector<XYZ> pts(point_count);
		for (int i = 0; i < point_count; ++i) pts[i] = ToCPP(points[i]);
		LN_NurbsCurve result;
		NurbsCurve::GlobalInterpolation(degree, pts, result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_global_interpolation_with_tangents(int degree, const XYZ_C* points, const XYZ_C* tangents,
		double tangent_factor, int point_count, LN_NurbsCurve_C* out_curve) {
		std::vector<XYZ> pts(point_count), tans(point_count);
		for (int i = 0; i < point_count; ++i) {
			pts[i] = ToCPP(points[i]);
			tans[i] = ToCPP(tangents[i]);
		}
		LN_NurbsCurve result;
		NurbsCurve::GlobalInterpolation(degree, pts, tans, tangent_factor, result);
		*out_curve = ToCAPI(result);
	}

	int nurbs_curve_cubic_local_interpolation(const XYZ_C* points, int point_count, LN_NurbsCurve_C* out_curve) {
		std::vector<XYZ> pts(point_count);
		for (int i = 0; i < point_count; ++i) pts[i] = ToCPP(points[i]);
		LN_NurbsCurve result;
		bool ok = NurbsCurve::CubicLocalInterpolation(pts, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	int nurbs_curve_least_squares_approximation(int degree, const XYZ_C* points, int point_count,
		int control_point_count, LN_NurbsCurve_C* out_curve) {
		std::vector<XYZ> pts(point_count);
		for (int i = 0; i < point_count; ++i) pts[i] = ToCPP(points[i]);
		LN_NurbsCurve result;
		bool ok = NurbsCurve::LeastSquaresApproximation(degree, pts, control_point_count, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	int nurbs_curve_weighted_constrained_least_squares(int degree, const XYZ_C* points, const double* point_weights,
		const XYZ_C* tangents, const int* tangent_indices, const double* tangent_weights,
		int tangent_count, int control_point_count, LN_NurbsCurve_C* out_curve) {
		std::vector<XYZ> pts, tans;
		std::vector<double> pws, tws;
		std::vector<int> tidx;
		for (int i = 0; i < control_point_count; ++i) {
			pts.push_back(ToCPP(points[i]));
			pws.push_back(point_weights[i]);
		}
		for (int i = 0; i < tangent_count; ++i) {
			tans.push_back(ToCPP(tangents[i]));
			tidx.push_back(tangent_indices[i]);
			tws.push_back(tangent_weights[i]);
		}
		LN_NurbsCurve result;
		bool ok = NurbsCurve::WeightedAndContrainedLeastSquaresApproximation(
			degree, pts, pws, tans, tidx, tws, control_point_count, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	void nurbs_curve_global_approximation_by_error_bound(int degree, const XYZ_C* points, int point_count,
		double max_error, LN_NurbsCurve_C* out_curve) {
		std::vector<XYZ> pts(point_count);
		for (int i = 0; i < point_count; ++i) pts[i] = ToCPP(points[i]);
		LN_NurbsCurve result;
		NurbsCurve::GlobalApproximationByErrorBound(degree, pts, max_error, result);
		*out_curve = ToCAPI(result);
	}

	XYZ_C nurbs_curve_get_point_on_curve(LN_NurbsCurve_C curve, double paramT) {
		return ToCAPI(NurbsCurve::GetPointOnCurve(FromCAPI(curve), paramT));
	}

	XYZ_C nurbs_curve_get_point_on_curve_by_corner_cut(LN_NurbsCurve_C curve, double paramT) {
		return ToCAPI(NurbsCurve::GetPointOnCurveByCornerCut(FromCAPI(curve), paramT));
	}

	int nurbs_curve_compute_rational_derivatives(LN_NurbsCurve_C curve, int derivative_order,
		double paramT, XYZ_C* out_derivatives) {
		auto derivatives = NurbsCurve::ComputeRationalCurveDerivatives(FromCAPI(curve), derivative_order, paramT);
		for (size_t i = 0; i < derivatives.size(); ++i) {
			out_derivatives[i] = ToCAPI(derivatives[i]);
		}
		return static_cast<int>(derivatives.size());
	}

	double nurbs_curve_curvature(LN_NurbsCurve_C curve, double paramT) {
		return NurbsCurve::Curvature(FromCAPI(curve), paramT);
	}

	double nurbs_curve_torsion(LN_NurbsCurve_C curve, double paramT) {
		return NurbsCurve::Torsion(FromCAPI(curve), paramT);
	}

	XYZ_C nurbs_curve_normal(LN_NurbsCurve_C curve, CurveNormal_C normal_type, double paramT) {
		return ToCAPI(NurbsCurve::Normal(FromCAPI(curve), static_cast<CurveNormal>(normal_type), paramT));
	}

	int nurbs_curve_project_normal(LN_NurbsCurve_C curve, XYZ_C* out_normals) {
		auto normals = NurbsCurve::ProjectNormal(FromCAPI(curve));
		for (size_t i = 0; i < normals.size(); ++i) {
			out_normals[i] = ToCAPI(normals[i]);
		}
		return static_cast<int>(normals.size());
	}

	double nurbs_curve_get_param_on_curve_by_point(LN_NurbsCurve_C curve, XYZ_C given_point) {
		return NurbsCurve::GetParamOnCurve(FromCAPI(curve), ToCPP(given_point));
	}

	double nurbs_curve_approximate_length(LN_NurbsCurve_C curve, IntegratorType_C integrator_type) {
		return NurbsCurve::ApproximateLength(FromCAPI(curve), static_cast<IntegratorType>(integrator_type));
	}

	double nurbs_curve_get_param_by_length(LN_NurbsCurve_C curve, double given_length, IntegratorType_C integrator_type) {
		return NurbsCurve::GetParamOnCurve(FromCAPI(curve), given_length, static_cast<IntegratorType>(integrator_type));
	}

	int nurbs_curve_get_params_by_equal_length(LN_NurbsCurve_C curve, double segment_length,
		IntegratorType_C integrator_type, double* out_params) {
		auto params = NurbsCurve::GetParamsOnCurve(FromCAPI(curve), segment_length,
			static_cast<IntegratorType>(integrator_type));
		for (size_t i = 0; i < params.size(); ++i) {
			out_params[i] = params[i];
		}
		return static_cast<int>(params.size());
	}

	int nurbs_curve_split_at(LN_NurbsCurve_C curve, double paramT,
		LN_NurbsCurve_C* out_left, LN_NurbsCurve_C* out_right) {
		LN_NurbsCurve left, right;
		bool ok = NurbsCurve::SplitAt(FromCAPI(curve), paramT, left, right);
		if (ok) {
			*out_left = ToCAPI(left);
			*out_right = ToCAPI(right);
		}
		return ok ? 1 : 0;
	}

	int nurbs_curve_segment(LN_NurbsCurve_C curve, double start_param, double end_param, LN_NurbsCurve_C* out_segment) {
		LN_NurbsCurve seg;
		bool ok = NurbsCurve::Segment(FromCAPI(curve), start_param, end_param, seg);
		if (ok) *out_segment = ToCAPI(seg);
		return ok ? 1 : 0;
	}

	int nurbs_curve_decompose_to_beziers(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_segments, int max_segments) {
		auto segments = NurbsCurve::DecomposeToBeziers(FromCAPI(curve));
		int n = std::min(static_cast<int>(segments.size()), max_segments);
		for (int i = 0; i < n; ++i) {
			out_segments[i] = ToCAPI(segments[i]);
		}
		return n;
	}

	int nurbs_curve_tessellate(LN_NurbsCurve_C curve, XYZ_C* out_points) {
		auto points = NurbsCurve::Tessellate(FromCAPI(curve));
		for (size_t i = 0; i < points.size(); ++i) {
			out_points[i] = ToCAPI(points[i]);
		}
		return static_cast<int>(points.size());
	}

	void nurbs_curve_create_transformed(LN_NurbsCurve_C curve, Matrix4d_C matrix, LN_NurbsCurve_C* out_curve) {
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

	void nurbs_curve_reverse(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::Reverse(FromCAPI(curve), result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_reparametrize_to_interval(LN_NurbsCurve_C curve, double min_val, double max_val, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::Reparametrize(FromCAPI(curve), min_val, max_val, result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_reparametrize_linear_rational(LN_NurbsCurve_C curve, double alpha, double beta,
		double gamma, double delta, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::Reparametrize(FromCAPI(curve), alpha, beta, gamma, delta, result);
		*out_curve = ToCAPI(result);
	}

	int nurbs_curve_insert_knot(LN_NurbsCurve_C curve, double knot_value, int times, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		int r = NurbsCurve::InsertKnot(FromCAPI(curve), knot_value, times, result);
		if (r >= 0) *out_curve = ToCAPI(result);
		return r;
	}

	int nurbs_curve_remove_knot(LN_NurbsCurve_C curve, double knot_value, int times, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		bool ok = NurbsCurve::RemoveKnot(FromCAPI(curve), knot_value, times, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	void nurbs_curve_remove_excessive_knots(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::RemoveExcessiveKnots(FromCAPI(curve), result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_refine_knot_vector(LN_NurbsCurve_C curve, const double* insert_knots, int insert_count, LN_NurbsCurve_C* out_curve) {
		std::vector<double> knots(insert_knots, insert_knots + insert_count);
		LN_NurbsCurve result;
		NurbsCurve::RefineKnotVector(FromCAPI(curve), knots, result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_elevate_degree(LN_NurbsCurve_C curve, int times, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::ElevateDegree(FromCAPI(curve), times, result);
		*out_curve = ToCAPI(result);
	}

	int nurbs_curve_reduce_degree(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		bool ok = NurbsCurve::ReduceDegree(FromCAPI(curve), result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	int nurbs_curve_is_closed(LN_NurbsCurve_C curve) {
		return NurbsCurve::IsClosed(FromCAPI(curve)) ? 1 : 0;
	}

	int nurbs_curve_is_linear(LN_NurbsCurve_C curve) {
		return NurbsCurve::IsLinear(FromCAPI(curve)) ? 1 : 0;
	}

	int nurbs_curve_is_clamped(LN_NurbsCurve_C curve) {
		return NurbsCurve::IsClamp(FromCAPI(curve)) ? 1 : 0;
	}

	int nurbs_curve_is_periodic(LN_NurbsCurve_C curve) {
		return NurbsCurve::IsPeriodic(FromCAPI(curve)) ? 1 : 0;
	}

	int nurbs_curve_can_compute_derivative(LN_NurbsCurve_C curve, double paramT) {
		return NurbsCurve::CanComputerDerivative(FromCAPI(curve), paramT) ? 1 : 0;
	}

	int nurbs_curve_control_point_reposition(LN_NurbsCurve_C curve, double paramT, int move_index,
		XYZ_C move_direction, double move_distance, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		bool ok = NurbsCurve::ControlPointReposition(FromCAPI(curve), paramT, move_index,
			ToCPP(move_direction), move_distance, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	void nurbs_curve_weight_modification(LN_NurbsCurve_C curve, double paramT, int move_index,
		double move_distance, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::WeightModification(FromCAPI(curve), paramT, move_index, move_distance, result);
		*out_curve = ToCAPI(result);
	}

	int nurbs_curve_neighbor_weights_modification(LN_NurbsCurve_C curve, double paramT, int move_index,
		double move_distance, double scale, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		bool ok = NurbsCurve::NeighborWeightsModification(FromCAPI(curve), paramT, move_index,
			move_distance, scale, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	void nurbs_curve_warping(LN_NurbsCurve_C curve, const double* warp_shape, int warp_shape_count,
		double warp_distance, XYZ_C plane_normal, double start_param, double end_param, LN_NurbsCurve_C* out_curve) {
		std::vector<double> shape(warp_shape, warp_shape + warp_shape_count);
		LN_NurbsCurve result;
		NurbsCurve::Warping(FromCAPI(curve), shape, warp_distance, ToCPP(plane_normal), start_param, end_param, result);
		*out_curve = ToCAPI(result);
	}

	int nurbs_curve_flattening(LN_NurbsCurve_C curve, XYZ_C line_start, XYZ_C line_end,
		double start_param, double end_param, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		bool ok = NurbsCurve::Flattening(FromCAPI(curve), ToCPP(line_start), ToCPP(line_end),
			start_param, end_param, result);
		if (ok) *out_curve = ToCAPI(result);
		return ok ? 1 : 0;
	}

	void nurbs_curve_bending(LN_NurbsCurve_C curve, double start_param, double end_param,
		XYZ_C bend_center, double radius, double cross_ratio, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::Bending(FromCAPI(curve), start_param, end_param, ToCPP(bend_center), radius, cross_ratio, result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_constraint_based_modification(LN_NurbsCurve_C curve,
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
		for (int i = 0; i < constraint_count; ++i) derivs[i] = ToCPP(derivative_constraints[i]);
		std::vector<int> aidx(applied_indices, applied_indices + constraint_count);
		std::vector<int> adegs(applied_degrees, applied_degrees + constraint_count);
		std::vector<int> fidx(fixed_cp_indices, fixed_cp_indices + fixed_count);
		LN_NurbsCurve result;
		NurbsCurve::ConstraintBasedModification(FromCAPI(curve), params, derivs, aidx, adegs, fidx, result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_to_clamp_curve(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::ToClampCurve(FromCAPI(curve), result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_to_unclamp_curve(LN_NurbsCurve_C curve, LN_NurbsCurve_C* out_curve) {
		LN_NurbsCurve result;
		NurbsCurve::ToUnclampCurve(FromCAPI(curve), result);
		*out_curve = ToCAPI(result);
	}

	void nurbs_curve_equally_tessellate(LN_NurbsCurve_C curve, XYZ_C* out_points, double* out_knots, int max_count) {
		std::vector<XYZ> pts;
		std::vector<double> knots;
		NurbsCurve::EquallyTessellate(FromCAPI(curve), pts, knots);
		int n = std::min(static_cast<int>(pts.size()), max_count);
		for (int i = 0; i < n; ++i) {
			out_points[i] = ToCAPI(pts[i]);
			out_knots[i] = knots[i];
		}
	}

}