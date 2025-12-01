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

    public static partial class LNLibAPI
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern LN_NurbsCurve nurbs_curve_create_line(XYZ start, XYZ end);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_create_arc(
            XYZ center, XYZ x_axis, XYZ y_axis,
            double start_rad, double end_rad,
            double x_radius, double y_radius,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_create_open_conic(
            XYZ start, XYZ start_tangent,
            XYZ end, XYZ end_tangent,
            XYZ point_on_conic,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_global_interpolation(
            int degree,
            IntPtr points,
            int point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_global_interpolation_with_tangents(
            int degree,
            IntPtr points,
            IntPtr tangents,
            double tangent_factor,
            int point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_cubic_local_interpolation(
            IntPtr points,
            int point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_least_squares_approximation(
            int degree,
            IntPtr points,
            int point_count,
            int control_point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_weighted_constrained_least_squares(
            int degree,
            IntPtr points,
            IntPtr point_weights,
            IntPtr tangents,
            IntPtr tangent_indices,
            IntPtr tangent_weights,
            int tangent_count,
            int control_point_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_global_approximation_by_error_bound(
            int degree,
            IntPtr points,
            int point_count,
            double max_error,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ nurbs_curve_get_point_on_curve(LN_NurbsCurve curve, double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ nurbs_curve_get_point_on_curve_by_corner_cut(LN_NurbsCurve curve, double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_compute_rational_derivatives(
            LN_NurbsCurve curve,
            int derivative_order,
            double paramT,
            IntPtr out_derivatives);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double nurbs_curve_curvature(LN_NurbsCurve curve, double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double nurbs_curve_torsion(LN_NurbsCurve curve, double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ nurbs_curve_normal(LN_NurbsCurve curve, CurveNormal normal_type, double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_project_normal(LN_NurbsCurve curve, IntPtr out_normals);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double nurbs_curve_get_param_on_curve_by_point(LN_NurbsCurve curve, XYZ given_point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double nurbs_curve_approximate_length(LN_NurbsCurve curve, IntegratorType integrator_type);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double nurbs_curve_get_param_by_length(LN_NurbsCurve curve, double given_length, IntegratorType integrator_type);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_get_params_by_equal_length(
            LN_NurbsCurve curve,
            double segment_length,
            IntegratorType integrator_type,
            IntPtr out_params);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_split_at(
            LN_NurbsCurve curve,
            double paramT,
            out LN_NurbsCurve out_left,
            out LN_NurbsCurve out_right);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_segment(
            LN_NurbsCurve curve,
            double start_param,
            double end_param,
            out LN_NurbsCurve out_segment);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_decompose_to_beziers(
            LN_NurbsCurve curve,
            IntPtr out_segments,
            int max_segments);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_tessellate(
            LN_NurbsCurve curve,
            IntPtr out_points);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_create_transformed(
            LN_NurbsCurve curve,
            Matrix4d matrix,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_reverse(LN_NurbsCurve curve, out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_reparametrize_to_interval(
            LN_NurbsCurve curve,
            double min_val,
            double max_val,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_reparametrize_linear_rational(
            LN_NurbsCurve curve,
            double alpha, double beta,
            double gamma, double delta,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_insert_knot(
            LN_NurbsCurve curve,
            double knot_value,
            int times,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_remove_knot(
            LN_NurbsCurve curve,
            double knot_value,
            int times,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_remove_excessive_knots(
            LN_NurbsCurve curve,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_refine_knot_vector(
            LN_NurbsCurve curve,
            IntPtr insert_knots,
            int insert_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_elevate_degree(
            LN_NurbsCurve curve,
            int times,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_reduce_degree(
            LN_NurbsCurve curve,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_is_closed(LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_is_linear(LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_is_clamped(LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_is_periodic(LN_NurbsCurve curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_can_compute_derivative(LN_NurbsCurve curve, double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_control_point_reposition(
            LN_NurbsCurve curve,
            double paramT,
            int move_index,
            XYZ move_direction,
            double move_distance,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_weight_modification(
            LN_NurbsCurve curve,
            double paramT,
            int move_index,
            double move_distance,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_neighbor_weights_modification(
            LN_NurbsCurve curve,
            double paramT,
            int move_index,
            double move_distance,
            double scale,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_warping(
            LN_NurbsCurve curve,
            IntPtr warp_shape,
            int warp_shape_count,
            double warp_distance,
            XYZ plane_normal,
            double start_param,
            double end_param,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int nurbs_curve_flattening(
            LN_NurbsCurve curve,
            XYZ line_start,
            XYZ line_end,
            double start_param,
            double end_param,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_bending(
            LN_NurbsCurve curve,
            double start_param,
            double end_param,
            XYZ bend_center,
            double radius,
            double cross_ratio,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_constraint_based_modification(
            LN_NurbsCurve curve,
            IntPtr constraint_params,
            IntPtr derivative_constraints,
            IntPtr applied_indices,
            IntPtr applied_degrees,
            IntPtr fixed_cp_indices,
            int constraint_count,
            int fixed_count,
            out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_to_clamp_curve(LN_NurbsCurve curve, out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_to_unclamp_curve(LN_NurbsCurve curve, out LN_NurbsCurve out_curve);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void nurbs_curve_equally_tessellate(
            LN_NurbsCurve curve,
            IntPtr out_points,
            IntPtr out_knots,
            int max_count);
    }
}