/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#pragma once
#include "LNLibDefinitions.h"
#include "UV_CAPI.h"
#include "XYZ_CAPI.h"
#include "XYZW_CAPI.h"
#include "LNObject_CAPI.h"
#include "NurbsCurve_CAPI.h"

#ifdef __cplusplus
extern "C" {
#endif

    typedef struct {
        int degree_u;
        int degree_v;
        const double* knot_vector_u;
        int knot_vector_count_u;
        const double* knot_vector_v;
        int knot_vector_count_v;
        const XYZW_C* control_points;
        int cp_rows;
        int cp_cols;
    } LN_NurbsSurface_C;

    LNLIB_EXPORT XYZ_C LNLIB_NURBSSUR_get_point_on_surface(LN_NurbsSurface_C surface, UV_C uv);
    LNLIB_EXPORT int LNLIB_NURBSSUR_compute_rational_surface_derivatives(
        LN_NurbsSurface_C surface,
        int derivative_order,
        UV_C uv,
        XYZ_C* out_derivatives);
    LNLIB_EXPORT void LNLIB_NURBSSUR_compute_rational_surface_first_order_derivative(
        LN_NurbsSurface_C surface,
        UV_C uv,
        XYZ_C* out_s,
        XYZ_C* out_su,
        XYZ_C* out_sv);
    LNLIB_EXPORT double LNLIB_NURBSSUR_curvature(LN_NurbsSurface_C surface, int curvature_type, UV_C uv);
    LNLIB_EXPORT XYZ_C LNLIB_NURBSSUR_Normal(LN_NurbsSurface_C surface, UV_C uv);

    LNLIB_EXPORT void LNLIB_NURBSSUR_swap(LN_NurbsSurface_C surface, LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_reverse(LN_NurbsSurface_C surface, int direction, LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_insert_knot(
        LN_NurbsSurface_C surface,
        double knot,
        int times,
        int is_u_direction,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_refine_knot_vector(
        LN_NurbsSurface_C surface,
        const double* insert_knots,
        int insert_count,
        int is_u_direction,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_decompose_to_beziers(
        LN_NurbsSurface_C surface,
        LN_NurbsSurface_C* out_patches,
        int max_patches);
    LNLIB_EXPORT void LNLIB_NURBSSUR_remove_knot(
        LN_NurbsSurface_C surface,
        double knot,
        int times,
        int is_u_direction,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_elevate_degree(
        LN_NurbsSurface_C surface,
        int times,
        int is_u_direction,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_reduce_degree(
        LN_NurbsSurface_C surface,
        int is_u_direction,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_equally_tessellate(
        LN_NurbsSurface_C surface,
        XYZ_C* out_points,
        UV_C* out_uvs,
        int max_count);
    LNLIB_EXPORT int LNLIB_NURBSSUR_is_closed(LN_NurbsSurface_C surface, int is_u_direction);
    LNLIB_EXPORT UV_C LNLIB_NURBSSUR_get_param_on_surface(LN_NurbsSurface_C surface, XYZ_C given_point);
    LNLIB_EXPORT UV_C LNLIB_NURBSSUR_get_param_on_surface_by_gsa(LN_NurbsSurface_C surface, XYZ_C given_point);
    LNLIB_EXPORT void LNLIB_NURBSSUR_reparametrize(
        LN_NurbsSurface_C surface,
        double min_u, double max_u,
        double min_v, double max_v,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_get_uv_tangent(
        LN_NurbsSurface_C surface,
        UV_C param,
        XYZ_C tangent,
        UV_C* out_uv_tangent);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_bilinear_surface(
        XYZ_C top_left, XYZ_C top_right,
        XYZ_C bottom_left, XYZ_C bottom_right,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_create_cylindrical_surface(
        XYZ_C origin, XYZ_C x_axis, XYZ_C y_axis,
        double start_rad, double end_rad,
        double radius, double height,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_ruled_surface(
        LN_NurbsCurve_C curve_0,
        LN_NurbsCurve_C curve_1,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_create_revolved_surface(
        XYZ_C origin, XYZ_C axis, double rad,
        LN_NurbsCurve_C profile,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_non_uniform_scaling_query_size(
        const XYZW_C* control_points,
        int rows,
        int cols,
        double x_factor,
        double y_factor,
        double z_factor,
        const XYZ_C reference_point,
        XYZW_C** output_control_points,
        int* out_rows,
        int* out_cols
    );
    LNLIB_EXPORT void LNLIB_NURBSSUR_make_corner_fillet_surface(const LN_NurbsCurve_C arc, LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_global_interpolation(
        const XYZ_C* points,
        int rows, int cols,
        int degree_u, int degree_v,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_bicubic_local_interpolation(
        const XYZ_C* points,
        int rows, int cols,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_global_approximation(
        const XYZ_C* points,
        int rows, int cols,
        int degree_u, int degree_v,
        int cp_rows, int cp_cols,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT int LNLIB_NURBSSUR_create_swung_surface(
        LN_NurbsCurve_C profile,
        LN_NurbsCurve_C trajectory,
        double scale,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_loft_surface(
        const LN_NurbsCurve_C* sections,
        int section_count,
        LN_NurbsSurface_C* out_surface,
        int custom_trajectory_degree,
        const double* custom_knot_vector,
        int knot_vector_count);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_generalized_translational_sweep_surface(
        LN_NurbsCurve_C profile,
        LN_NurbsCurve_C trajectory,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_sweep_surface(
        LN_NurbsCurve_C profile,
        LN_NurbsCurve_C trajectory,
        int min_profiles,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_sweep_surface_by_trajectory_degree(
        LN_NurbsCurve_C profile,
        LN_NurbsCurve_C trajectory,
        int min_profiles,
        int trajectory_degree,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_gordon_surface(
        const LN_NurbsCurve_C* u_curves,
        int u_count,
        const LN_NurbsCurve_C* v_curves,
        int v_count,
        const XYZ_C* intersections,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT void LNLIB_NURBSSUR_create_coons_surface(
        LN_NurbsCurve_C left,
        LN_NurbsCurve_C bottom,
        LN_NurbsCurve_C right,
        LN_NurbsCurve_C top,
        LN_NurbsSurface_C* out_surface);
    LNLIB_EXPORT double LNLIB_NURBSSUR_approximate_area(LN_NurbsSurface_C surface, int integrator_type);
    LNLIB_EXPORT LNLIB_Mesh_C LNLIB_NURBSSUR_triangulate(
        LN_NurbsSurface_C surface,
        int resolution_u,
        int resolution_v,
        int use_delaunay);

#ifdef __cplusplus
}
#endif