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
    public struct LN_NurbsSurface
    {
        public int degree_u;
        public int degree_v;
        public IntPtr knot_vector_u;
        public int knot_vector_count_u;
        public IntPtr knot_vector_v;
        public int knot_vector_count_v;
        public IntPtr control_points;
        public int cp_rows;
        public int cp_cols;
    }

    public static partial class LNLibNurbsSurface
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_get_point_on_surface")]
        public static extern XYZ LNLIB_NURBSSUR_get_point_on_surface(
        LN_NurbsSurface surface,
        UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_compute_rational_surface_derivatives")]
        public static extern int LNLIB_NURBSSUR_compute_rational_surface_derivatives(
            LN_NurbsSurface surface,
            int derivative_order,
            UV uv,
            [Out] XYZ[] out_derivatives);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_compute_rational_surface_first_order_derivative")]
        public static extern void LNLIB_NURBSSUR_compute_rational_surface_first_order_derivative(
            LN_NurbsSurface surface,
            UV uv,
            out XYZ out_s,
            out XYZ out_su,
            out XYZ out_sv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_curvature")]
        public static extern double LNLIB_NURBSSUR_curvature(
            LN_NurbsSurface surface,
            SurfaceCurvature curvature_type,
            UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_Normal")]
        public static extern XYZ LNLIB_NURBSSUR_Normal(
            LN_NurbsSurface surface,
            UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_swap")]
        public static extern void LNLIB_NURBSSUR_swap(
            LN_NurbsSurface surface,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_reverse")]
        public static extern void LNLIB_NURBSSUR_reverse(
            LN_NurbsSurface surface,
            SurfaceDirection direction,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_insert_knot")]
        public static extern void LNLIB_NURBSSUR_insert_knot(
            LN_NurbsSurface surface,
            double knot,
            int times,
            int is_u_direction,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_refine_knot_vector")]
        public static extern void LNLIB_NURBSSUR_refine_knot_vector(
            LN_NurbsSurface surface,
            [In] double[] insert_knots,
            int insert_count,
            int is_u_direction,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_decompose_to_beziers")]
        public static extern int LNLIB_NURBSSUR_decompose_to_beziers(
            LN_NurbsSurface surface,
            [Out] LN_NurbsSurface[] out_patches,
            int max_patches);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_remove_knot")]
        public static extern void LNLIB_NURBSSUR_remove_knot(
            LN_NurbsSurface surface,
            double knot,
            int times,
            int is_u_direction,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_elevate_degree")]
        public static extern void LNLIB_NURBSSUR_elevate_degree(
            LN_NurbsSurface surface,
            int times,
            int is_u_direction,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_reduce_degree")]
        public static extern int LNLIB_NURBSSUR_reduce_degree(
            LN_NurbsSurface surface,
            int is_u_direction,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_equally_tessellate")]
        public static extern void LNLIB_NURBSSUR_equally_tessellate(
            LN_NurbsSurface surface,
            [Out] XYZ[] out_points,
            [Out] UV[] out_uvs,
            int max_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_is_closed")]
        public static extern int LNLIB_NURBSSUR_is_closed(
            LN_NurbsSurface surface,
            int is_u_direction);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_get_param_on_surface")]
        public static extern UV LNLIB_NURBSSUR_get_param_on_surface(
            LN_NurbsSurface surface,
            XYZ given_point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_get_param_on_surface_by_gsa")]
        public static extern UV LNLIB_NURBSSUR_get_param_on_surface_by_gsa(
            LN_NurbsSurface surface,
            XYZ given_point);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_reparametrize")]
        public static extern void LNLIB_NURBSSUR_reparametrize(
            LN_NurbsSurface surface,
            double min_u,
            double max_u,
            double min_v,
            double max_v,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_get_uv_tangent")]
        public static extern int LNLIB_NURBSSUR_get_uv_tangent(
            LN_NurbsSurface surface,
            UV param,
            XYZ tangent,
            out UV out_uv_tangent);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_bilinear_surface")]
        public static extern void LNLIB_NURBSSUR_create_bilinear_surface(
            XYZ top_left,
            XYZ top_right,
            XYZ bottom_left,
            XYZ bottom_right,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_cylindrical_surface")]
        public static extern int LNLIB_NURBSSUR_create_cylindrical_surface(
            XYZ origin,
            XYZ x_axis,
            XYZ y_axis,
            double start_rad,
            double end_rad,
            double radius,
            double height,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_ruled_surface")]
        public static extern void LNLIB_NURBSSUR_create_ruled_surface(
            LN_NurbsCurve curve_0,
            LN_NurbsCurve curve_1,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_revolved_surface")]
        public static extern int LNLIB_NURBSSUR_create_revolved_surface(
            XYZ origin,
            XYZ axis,
            double rad,
            LN_NurbsCurve profile,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_non_uniform_scaling_query_size")]
        public static extern void LNLIB_NURBSSUR_non_uniform_scaling_query_size(
            [In] XYZW[] control_points,
            int rows,
            int cols,
            double x_factor,
            double y_factor,
            double z_factor,
            XYZ reference_point,
            out IntPtr output_control_points,
            out int out_rows,
            out int out_cols);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_make_corner_fillet_surface")]
        public static extern void LNLIB_NURBSSUR_make_corner_fillet_surface(
            LN_NurbsCurve arc,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_global_interpolation")]
        public static extern void LNLIB_NURBSSUR_global_interpolation(
            [In] XYZ[] points,
            int rows,
            int cols,
            int degree_u,
            int degree_v,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_bicubic_local_interpolation")]
        public static extern int LNLIB_NURBSSUR_bicubic_local_interpolation(
            [In] XYZ[] points,
            int rows,
            int cols,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_global_approximation")]
        public static extern int LNLIB_NURBSSUR_global_approximation(
            [In] XYZ[] points,
            int rows,
            int cols,
            int degree_u,
            int degree_v,
            int cp_rows,
            int cp_cols,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_swung_surface")]
        public static extern int LNLIB_NURBSSUR_create_swung_surface(
            LN_NurbsCurve profile,
            LN_NurbsCurve trajectory,
            double scale,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_loft_surface")]
        public static extern void LNLIB_NURBSSUR_create_loft_surface(
            [In] LN_NurbsCurve[] sections,
            int section_count,
            out LN_NurbsSurface out_surface,
            int custom_trajectory_degree,
            [In] double[] custom_knot_vector,
            int knot_vector_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_generalized_translational_sweep_surface")]
        public static extern void LNLIB_NURBSSUR_create_generalized_translational_sweep_surface(
            LN_NurbsCurve profile,
            LN_NurbsCurve trajectory,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_sweep_surface")]
        public static extern void LNLIB_NURBSSUR_create_sweep_surface(
            LN_NurbsCurve profile,
            LN_NurbsCurve trajectory,
            int min_profiles,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_sweep_surface_by_trajectory_degree")]
        public static extern void LNLIB_NURBSSUR_create_sweep_surface_by_trajectory_degree(
            LN_NurbsCurve profile,
            LN_NurbsCurve trajectory,
            int min_profiles,
            int trajectory_degree,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_gordon_surface")]
        public static extern void LNLIB_NURBSSUR_create_gordon_surface(
            [In] LN_NurbsCurve[] u_curves,
            int u_count,
            [In] LN_NurbsCurve[] v_curves,
            int v_count,
            [In] XYZ[] intersections,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_create_coons_surface")]
        public static extern void LNLIB_NURBSSUR_create_coons_surface(
            LN_NurbsCurve left,
            LN_NurbsCurve bottom,
            LN_NurbsCurve right,
            LN_NurbsCurve top,
            out LN_NurbsSurface out_surface);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_approximate_area")]
        public static extern double LNLIB_NURBSSUR_approximate_area(
            LN_NurbsSurface surface,
            IntegratorType integrator_type);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_NURBSSUR_triangulate")]
        public static extern LN_Mesh LNLIB_NURBSSUR_triangulate(
            LN_NurbsSurface surface,
            int resolution_u,
            int resolution_v,
            int use_delaunay);
    }
}