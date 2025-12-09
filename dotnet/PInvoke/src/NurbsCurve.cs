/*
 * Author:
 * 2025/11/30 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System;
using System.Runtime.InteropServices;
using static LNLibSharp.LNLibDefinitions;

namespace LNLibSharp
{
    [StructLayout(LayoutKind.Sequential)]
    public struct LN_NurbsCurve
    {
        public int degree;
        public IntPtr knot_vector;
        public int knot_count;
        public IntPtr control_points;
        public int control_point_count;
    }

    public static partial class LNLibNurbsCurve
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_get_point_on_curve")]
        public static extern XYZ GetPointOnCurve(
        LN_NurbsCurve curve,
        double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_compute_rational_curve_derivatives")]
        public static extern int ComputeRationalCurveDerivatives(
            LN_NurbsCurve curve,
            int derivative_order,
            double param_t,
            [Out] XYZ[] out_derivatives);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_can_compute_derivative")]
        public static extern int CanComputeDerivative(
            LN_NurbsCurve curve,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_curvature")]
        public static extern double Curvature(
            LN_NurbsCurve curve,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_torsion")]
        public static extern double Torsion(
            LN_NurbsCurve curve,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_insert_knot")]
        public static extern int InsertKnot(
            LN_NurbsCurve curve,
            double knot,
            int times,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_get_point_on_curve_by_corner_cut")]
        public static extern XYZ GetPointOnCurveByCornerCut(
            LN_NurbsCurve curve,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_refine_knot_vector")]
        public static extern void RefineKnotVector(
            LN_NurbsCurve curve,
            [In] double[] insert_knots,
            int insert_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_decompose_to_beziers")]
        public static extern int DecomposeToBeziers(
            LN_NurbsCurve curve,
            [Out] LN_NurbsCurve[] out_segments,
            int max_segments);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_remove_knot")]
        public static extern int RemoveKnot(
            LN_NurbsCurve curve,
            double knot,
            int times,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_remove_excessive_knots")]
        public static extern void RemoveExcessiveKnots(
            LN_NurbsCurve curve,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_elevate_degree")]
        public static extern void ElevateDegree(
            LN_NurbsCurve curve,
            int times,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_reduce_degree")]
        public static extern int ReduceDegree(
            LN_NurbsCurve curve,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_equally_tessellate")]
        public static extern void EquallyTessellate(
            LN_NurbsCurve curve,
            int max_count,
            [Out] XYZ[] out_points,
            [Out] double[] out_knots);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_is_closed")]
        public static extern int IsClosed(
            LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_get_param_on_curve")]
        public static extern double GetParamOnCurve(
            LN_NurbsCurve curve,
            XYZ given_point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_create_transformed")]
        public static extern void CreateTransformed(
            LN_NurbsCurve curve,
            Matrix4d matrix,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_reparametrize")]
        public static extern void Reparametrize(
            LN_NurbsCurve curve,
            double min_val,
            double max_val,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_reparametrize_linear_rational")]
        public static extern void ReparametrizeLinearRational(
            LN_NurbsCurve curve,
            double alpha,
            double beta,
            double gamma,
            double delta,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_reverse")]
        public static extern void Reverse(
            LN_NurbsCurve curve,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_split_at")]
        public static extern int SplitAt(
            LN_NurbsCurve curve,
            double param_t,
            out LN_NurbsCurve out_left,
            out LN_NurbsCurve out_right);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_segment")]
        public static extern int Segment(
            LN_NurbsCurve curve,
            double start_param,
            double end_param,
            out LN_NurbsCurve out_segment);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_merge")]
        public static extern int Merge(
            LN_NurbsCurve left,
            LN_NurbsCurve right,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_offset")]
        public static extern void Offset(
            LN_NurbsCurve curve,
            double offset,
            OffsetType type,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_create_line")]
        public static extern LN_NurbsCurve CreateLine(
            XYZ start,
            XYZ end);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_create_cubic_hermite")]
        public static extern void CreateCubicHermite(
            [In] XYZ[] through_points,
            int through_points_count,
            [In] XYZ[] tangents,
            int tangents_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_create_arc")]
        public static extern int CreateArc(
            XYZ center,
            XYZ x_axis,
            XYZ y_axis,
            double start_rad,
            double end_rad,
            double x_radius,
            double y_radius,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_create_one_conic_arc")]
        public static extern bool CreateOneConicArc(
            XYZ start,
            XYZ start_tangent,
            XYZ end,
            XYZ end_tangent,
            XYZ point_on_conic,
            out XYZ out_project_point,
            out double out_project_point_weight);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_split_arc")]
        public static extern void SplitArc(
            XYZ start,
            XYZ project_point,
            double project_point_weight,
            XYZ end,
            out XYZ out_insert_point_at_start_side,
            out XYZ out_split_point,
            out XYZ out_insert_point_at_end_side,
            out double out_insert_weight);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_create_open_conic")]
        public static extern int CreateOpenConic(
            XYZ start,
            XYZ start_tangent,
            XYZ end,
            XYZ end_tangent,
            XYZ point_on_conic,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_global_interpolation")]
        public static extern void GlobalInterpolation(
            int degree,
            [In] XYZ[] points,
            int point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_global_interpolation_with_tangents")]
        public static extern void GlobalInterpolationWthTangents(
            int degree,
            [In] XYZ[] points,
            [In] XYZ[] tangents,
            double tangent_factor,
            int point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_cubic_local_interpolation")]
        public static extern int CubicLocalInterpolation(
            [In] XYZ[] points,
            int point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_least_squares_approximation")]
        public static extern int LeastSquaresApproximation(
            int degree,
            [In] XYZ[] points,
            int point_count,
            int control_point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_weighted_constrained_least_squares")]
        public static extern int WeightedConstrainedLeastSquares(
            int degree,
            [In] XYZ[] points,
            [In] double[] point_weights,
            [In] XYZ[] tangents,
            [In] int[] tangent_indices,
            [In] double[] tangent_weights,
            int tangent_count,
            int control_point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_computer_remove_knot_error_bound")]
        public static extern double ComputerRemoveKnotErrorBound(
            LN_NurbsCurve curve,
            int removal_index);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_remove_knots_by_given_bound")]
        public static extern void RemoveKnotsByGivenBound(
            LN_NurbsCurve curve,
            [In] double[] parameters,
            int params_count,
            double maxError,
            [Out] double[] out_errors,
            [In, Out] ref int out_errors_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_global_approximation_by_error_bound")]
        public static extern void GlobalApproximationByErrorBound(
            int degree,
            [In] XYZ[] points,
            int point_count,
            double max_error,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_fit_with_conic")]
        public static extern bool FitWithConic(
            [In] XYZ[] through_points,
            int through_points_count,
            int start_point_index,
            int end_point_index,
            XYZ start_tangent,
            XYZ end_tangent,
            double max_error,
            [Out] XYZW[] middle_control_points,
            [In, Out] ref int out_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_fit_with_cubic")]
        public static extern int FitWithCubic(
            [In] XYZ[] through_points,
            int through_points_count,
            int start_point_index,
            int end_point_index,
            XYZ start_tangent,
            XYZ end_tangent,
            double max_error,
            [Out] XYZW[] middle_control_points,
            [In, Out] ref int out_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_normal")]
        public static extern XYZ Normal(
            LN_NurbsCurve curve,
            CurveNormal normal_type,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_project_normal")]
        public static extern int ProjectNormal(
            LN_NurbsCurve curve,
            [Out] XYZ[] out_normals);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_control_point_reposition")]
        public static extern int ControlPointReposition(
            LN_NurbsCurve curve,
            double param_t,
            int move_index,
            XYZ move_direction,
            double move_distance,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_weight_modification")]
        public static extern void WeightModification(
            LN_NurbsCurve curve,
            double param_t,
            int move_index,
            double move_distance,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_neighbor_weights_modification")]
        public static extern int NeighborWeightsModification(
            LN_NurbsCurve curve,
            double param_t,
            int move_index,
            double move_distance,
            double scale,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_warping")]
        public static extern void Warping(
            LN_NurbsCurve curve,
            [In] double[] warp_shape,
            int warp_shape_count,
            double warp_distance,
            XYZ plane_normal,
            double start_param,
            double end_param,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_flattening")]
        public static extern int Flattening(
            LN_NurbsCurve curve,
            XYZ line_start,
            XYZ line_end,
            double start_param,
            double end_param,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_bending")]
        public static extern void Bending(
            LN_NurbsCurve curve,
            double start_param,
            double end_param,
            XYZ bend_center,
            double radius,
            double cross_ratio,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_constraint_based_modification")]
        public static extern void ConstraintBasedModification(
            LN_NurbsCurve curve,
            [In] double[] constraint_params,
            [In] XYZ[] derivative_constraints,
            [In] int[] applied_indices,
            [In] int[] applied_degrees,
            [In] int[] fixed_cp_indices,
            int constraint_count,
            int fixed_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_is_clamp")]
        public static extern int IsClamp(
            LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_to_clamp_curve")]
        public static extern void ToClampCurve(
            LN_NurbsCurve curve,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_is_periodic")]
        public static extern int IsPeriodic(
            LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_to_unclamp_curve")]
        public static extern void ToUnclampCurve(
            LN_NurbsCurve curve,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_is_linear")]
        public static extern int IsLinear(
            LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_is_arc")]
        public static extern int IsArc(
            LN_NurbsCurve curve,
            IntPtr out_arc_info);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_approximate_length")]
        public static extern double ApproximateLength(
            LN_NurbsCurve curve,
            IntegratorType integrator_type);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_get_param_by_length")]
        public static extern double GetParamByLength(
            LN_NurbsCurve curve,
            double given_length,
            IntegratorType integrator_type);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_get_params_by_equal_length")]
        public static extern int GetParamsByEqualLength(
            LN_NurbsCurve curve,
            double segment_length,
            IntegratorType integrator_type,
            [Out] double[] out_params);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSCUR_tessellate")]
        public static extern int Tessellate(
            LN_NurbsCurve curve,
            [Out] XYZ[] out_points);
    }
}