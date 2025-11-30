/*
 * Author:
 * 2025/11/30 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System.Runtime.InteropServices;

namespace LNLibSharp
{
[StructLayout(LayoutKind.Sequential)]
public struct LN_NurbsSurface
{
    public int degree_u;
    public IntPtr knot_vector_u; // double[]
    public int knot_count_u;
    public int degree_v;
    public IntPtr knot_vector_v; // double[]
    public int knot_count_v;
    public IntPtr control_points; // XYZW[]
    public int control_point_rows;
    public int control_point_cols;
}

public static partial class LNLibAPI
{
    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZ nurbs_surface_get_point_on_surface(LN_NurbsSurface surface, UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_compute_rational_derivatives(
        LN_NurbsSurface surface,
        int derivative_order,
        UV uv,
        IntPtr out_derivatives);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_compute_first_order_derivative(
        LN_NurbsSurface surface,
        UV uv,
        out XYZ out_S,
        out XYZ out_Su,
        out XYZ out_Sv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double nurbs_surface_curvature(LN_NurbsSurface surface, int curvature_type, UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZ nurbs_surface_normal(LN_NurbsSurface surface, UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_swap_uv(LN_NurbsSurface surface, out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_reverse(LN_NurbsSurface surface, int direction, out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_is_closed(LN_NurbsSurface surface, int is_u_direction);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_insert_knot(
        LN_NurbsSurface surface,
        double knot_value,
        int times,
        int is_u_direction,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_refine_knot_vector(
        LN_NurbsSurface surface,
        IntPtr insert_knots,
        int insert_count,
        int is_u_direction,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_remove_knot(
        LN_NurbsSurface surface,
        double knot_value,
        int times,
        int is_u_direction,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_elevate_degree(
        LN_NurbsSurface surface,
        int times,
        int is_u_direction,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_reduce_degree(
        LN_NurbsSurface surface,
        int is_u_direction,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_decompose_to_beziers(
        LN_NurbsSurface surface,
        IntPtr out_patches,
        int max_patches);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_equally_tessellate(
        LN_NurbsSurface surface,
        IntPtr out_points,
        IntPtr out_uvs,
        int max_count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV nurbs_surface_get_param_on_surface(LN_NurbsSurface surface, XYZ given_point);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV nurbs_surface_get_param_on_surface_by_gsa(LN_NurbsSurface surface, XYZ given_point);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_get_uv_tangent(
        LN_NurbsSurface surface,
        UV param,
        XYZ tangent,
        out UV out_uv_tangent);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_reparametrize(
        LN_NurbsSurface surface,
        double min_u, double max_u,
        double min_v, double max_v,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_bilinear(
        XYZ top_left, XYZ top_right,
        XYZ bottom_left, XYZ bottom_right,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_create_cylindrical(
        XYZ origin, XYZ x_axis, XYZ y_axis,
        double start_rad, double end_rad,
        double radius, double height,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_ruled(
        LN_NurbsCurve curve0,
        LN_NurbsCurve curve1,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_create_revolved(
        XYZ origin, XYZ axis, double rad,
        LN_NurbsCurve profile,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_global_interpolation(
        IntPtr points,
        int rows, int cols,
        int degree_u, int degree_v,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_bicubic_local_interpolation(
        IntPtr points,
        int rows, int cols,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_global_approximation(
        IntPtr points,
        int rows, int cols,
        int degree_u, int degree_v,
        int ctrl_rows, int ctrl_cols,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int nurbs_surface_create_swung(
        LN_NurbsCurve profile,
        LN_NurbsCurve trajectory,
        double scale,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_loft(
        IntPtr sections,
        int section_count,
        out LN_NurbsSurface out_surface,
        int custom_trajectory_degree,
        IntPtr custom_knots,
        int knot_count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_generalized_translational_sweep(
        LN_NurbsCurve profile,
        LN_NurbsCurve trajectory,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_sweep_interpolated(
        LN_NurbsCurve profile,
        LN_NurbsCurve trajectory,
        int min_profiles,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_sweep_noninterpolated(
        LN_NurbsCurve profile,
        LN_NurbsCurve trajectory,
        int min_profiles,
        int trajectory_degree,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_gordon(
        IntPtr u_curves,
        int u_count,
        IntPtr v_curves,
        int v_count,
        IntPtr intersections,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void nurbs_surface_create_coons(
        LN_NurbsCurve left,
        LN_NurbsCurve bottom,
        LN_NurbsCurve right,
        LN_NurbsCurve top,
        out LN_NurbsSurface out_surface);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double nurbs_surface_approximate_area(LN_NurbsSurface surface, int integrator_type);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern LN_Mesh nurbs_surface_triangulate(
        LN_NurbsSurface surface,
        int resolution_u,
        int resolution_v,
        int use_delaunay);
}
}