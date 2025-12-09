/*
 * Author:
 * 2025/11/27 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 *
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

#include "NurbsSurface_CAPI.h"
#include "NurbsSurface.h"
#include "NurbsCurve_CAPI.h"
#include "LNLibDefinitions.h"
#include <vector>
#include <thread>

using namespace LNLib;

static thread_local std::vector<double> g_knot_u_buffer;
static thread_local std::vector<double> g_knot_v_buffer;
static thread_local std::vector<std::vector<XYZW>> g_cp_2d_buffer;
static thread_local std::vector<XYZW> g_cp_flat_buffer;
static thread_local std::vector<LN_NurbsSurface> g_surface_buffer;
static thread_local std::vector<XYZ> g_xyz_buffer;
static thread_local std::vector<UV> g_uv_buffer;

static LN_NurbsSurface FromCAPI(LN_NurbsSurface_C c) {
	LN_NurbsSurface s;
	s.DegreeU = c.degree_u;
	s.DegreeV = c.degree_v;
	s.KnotVectorU.assign(c.knot_vector_u, c.knot_vector_u + c.knot_vector_count_u);
	s.KnotVectorV.assign(c.knot_vector_v, c.knot_vector_v + c.knot_vector_count_v);
	s.ControlPoints.resize(c.cp_rows, std::vector<XYZW>(c.cp_cols));
	for (int i = 0; i < c.cp_rows; ++i) {
		for (int j = 0; j < c.cp_cols; ++j) {
			XYZW_C cp = c.control_points[i * c.cp_cols + j];
			s.ControlPoints[i][j] = XYZW(cp.wx, cp.wy, cp.wz, cp.w);
		}
	}
	return s;
}

static LN_NurbsSurface_C ToCAPI(const LN_NurbsSurface& s) {
	g_knot_u_buffer = s.KnotVectorU;
	g_knot_v_buffer = s.KnotVectorV;
	g_cp_2d_buffer = s.ControlPoints;
	int rows = static_cast<int>(g_cp_2d_buffer.size());
	int cols = rows > 0 ? static_cast<int>(g_cp_2d_buffer[0].size()) : 0;
	g_cp_flat_buffer.clear();
	g_cp_flat_buffer.reserve(rows * cols);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			g_cp_flat_buffer.push_back(g_cp_2d_buffer[i][j]);
		}
	}

	LN_NurbsSurface_C c;
	c.degree_u = s.DegreeU;
	c.degree_v = s.DegreeV;
	c.knot_vector_u = g_knot_u_buffer.data();
	c.knot_vector_count_u = static_cast<int>(g_knot_u_buffer.size());
	c.knot_vector_v = g_knot_v_buffer.data();
	c.knot_vector_count_v = static_cast<int>(g_knot_v_buffer.size());
	c.control_points = reinterpret_cast<const XYZW_C*>(g_cp_flat_buffer.data());
	c.cp_rows = rows;
	c.cp_cols = cols;
	return c;
}

XYZ_C LNLIB_NURBSSUR_get_point_on_surface(LN_NurbsSurface_C surface, UV_C uv) {
	return FromXYZ(NurbsSurface::GetPointOnSurface(FromCAPI(surface), ToUV(uv)));
}

int LNLIB_NURBSSUR_compute_rational_surface_derivatives(LN_NurbsSurface_C surface, int derivative_order,
	UV_C uv, XYZ_C* out_derivatives) {
	auto result = NurbsSurface::ComputeRationalSurfaceDerivatives(FromCAPI(surface), derivative_order, ToUV(uv));
	int count = 0;
	for (size_t i = 0; i < result.size(); ++i) {
		for (size_t j = 0; j < result[i].size(); ++j) {
			out_derivatives[count++] = FromXYZ(result[i][j]);
		}
	}
	return count;
}

void LNLIB_NURBSSUR_compute_rational_surface_first_order_derivative(LN_NurbsSurface_C surface, UV_C uv,
	XYZ_C* out_S, XYZ_C* out_Su, XYZ_C* out_Sv) {
	XYZ S, Su, Sv;
	NurbsSurface::ComputeRationalSurfaceFirstOrderDerivative(FromCAPI(surface), ToUV(uv), S, Su, Sv);
	*out_S = FromXYZ(S);
	*out_Su = FromXYZ(Su);
	*out_Sv = FromXYZ(Sv);
}

double LNLIB_NURBSSUR_curvature(LN_NurbsSurface_C surface, int curvature_type, UV_C uv) {
	return NurbsSurface::Curvature(FromCAPI(surface), static_cast<SurfaceCurvature>(curvature_type), ToUV(uv));
}

XYZ_C LNLIB_NURBSSUR_Normal(LN_NurbsSurface_C surface, UV_C uv) {
	return FromXYZ(NurbsSurface::Normal(FromCAPI(surface), ToUV(uv)));
}

void LNLIB_NURBSSUR_swap(LN_NurbsSurface_C surface, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::Swap(FromCAPI(surface), result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_reverse(LN_NurbsSurface_C surface, int direction, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::Reverse(FromCAPI(surface), static_cast<SurfaceDirection>(direction), result);
	*out_surface = ToCAPI(result);
}

int LNLIB_NURBSSUR_is_closed(LN_NurbsSurface_C surface, int is_u_direction) {
	return NurbsSurface::IsClosed(FromCAPI(surface), is_u_direction != 0) ? 1 : 0;
}

void LNLIB_NURBSSUR_insert_knot(LN_NurbsSurface_C surface, double knot_value, int times,
	int is_u_direction, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::InsertKnot(FromCAPI(surface), knot_value, times, is_u_direction != 0, result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_refine_knot_vector(LN_NurbsSurface_C surface, const double* insert_knots,
	int insert_count, int is_u_direction, LN_NurbsSurface_C* out_surface) {
	std::vector<double> knots(insert_knots, insert_knots + insert_count);
	LN_NurbsSurface result;
	NurbsSurface::RefineKnotVector(FromCAPI(surface), knots, is_u_direction != 0, result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_remove_knot(LN_NurbsSurface_C surface, double knot_value, int times,
	int is_u_direction, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::RemoveKnot(FromCAPI(surface), knot_value, times, is_u_direction != 0, result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_elevate_degree(LN_NurbsSurface_C surface, int times, int is_u_direction,
	LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::ElevateDegree(FromCAPI(surface), times, is_u_direction != 0, result);
	*out_surface = ToCAPI(result);
}

int LNLIB_NURBSSUR_reduce_degree(LN_NurbsSurface_C surface, int is_u_direction, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	bool ok = NurbsSurface::ReduceDegree(FromCAPI(surface), is_u_direction != 0, result);
	if (ok) *out_surface = ToCAPI(result);
	return ok ? 1 : 0;
}

int LNLIB_NURBSSUR_decompose_to_beziers(LN_NurbsSurface_C surface, LN_NurbsSurface_C* out_patches, int max_patches) {
	auto patches = NurbsSurface::DecomposeToBeziers(FromCAPI(surface));
	int n = std::min(static_cast<int>(patches.size()), max_patches);
	for (int i = 0; i < n; ++i) {
		out_patches[i] = ToCAPI(patches[i]);
	}
	return n;
}

void LNLIB_NURBSSUR_equally_tessellate(LN_NurbsSurface_C surface, XYZ_C* out_points, UV_C* out_uvs, int max_count) {
	std::vector<XYZ> pts;
	std::vector<UV> uvs;
	NurbsSurface::EquallyTessellate(FromCAPI(surface), pts, uvs);
	int n = std::min(static_cast<int>(pts.size()), max_count);
	for (int i = 0; i < n; ++i) {
		out_points[i] = FromXYZ(pts[i]);
		out_uvs[i] = FromUV(uvs[i]);
	}
}

UV_C LNLIB_NURBSSUR_get_param_on_surface(LN_NurbsSurface_C surface, XYZ_C given_point) {
	return FromUV(NurbsSurface::GetParamOnSurface(FromCAPI(surface), ToXYZ(given_point)));
}

UV_C LNLIB_NURBSSUR_get_param_on_surface_by_gsa(LN_NurbsSurface_C surface, XYZ_C given_point) {
	return FromUV(NurbsSurface::GetParamOnSurfaceByGSA(FromCAPI(surface), ToXYZ(given_point)));
}

int LNLIB_NURBSSUR_get_uv_tangent(LN_NurbsSurface_C surface, UV_C param, XYZ_C tangent, UV_C* out_uv_tangent) {
	UV uvTangent;
	bool ok = NurbsSurface::GetUVTangent(FromCAPI(surface), ToUV(param), ToXYZ(tangent), uvTangent);
	if (ok) *out_uv_tangent = FromUV(uvTangent);
	return ok ? 1 : 0;
}

void LNLIB_NURBSSUR_reparametrize(LN_NurbsSurface_C surface, double min_u, double max_u,
	double min_v, double max_v, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::Reparametrize(FromCAPI(surface), min_u, max_u, min_v, max_v, result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_create_bilinear_surface(XYZ_C top_left, XYZ_C top_right, XYZ_C bottom_left, XYZ_C bottom_right,
	LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::CreateBilinearSurface(ToXYZ(top_left), ToXYZ(top_right), ToXYZ(bottom_left), ToXYZ(bottom_right), result);
	*out_surface = ToCAPI(result);
}

int LNLIB_NURBSSUR_create_cylindrical_surface(XYZ_C origin, XYZ_C x_axis, XYZ_C y_axis,
	double start_rad, double end_rad, double radius, double height, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	bool ok = NurbsSurface::CreateCylindricalSurface(ToXYZ(origin), ToXYZ(x_axis), ToXYZ(y_axis),
		start_rad, end_rad, radius, height, result);
	if (ok) *out_surface = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSSUR_create_ruled_surface(LN_NurbsCurve_C curve0, LN_NurbsCurve_C curve1, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::CreateRuledSurface(FromCAPI(curve0), FromCAPI(curve1), result);
	*out_surface = ToCAPI(result);
}

int LNLIB_NURBSSUR_create_revolved_surface(XYZ_C origin, XYZ_C axis, double rad, LN_NurbsCurve_C profile,
	LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	bool ok = NurbsSurface::CreateRevolvedSurface(ToXYZ(origin), ToXYZ(axis), rad, FromCAPI(profile), result);
	if (ok) *out_surface = ToCAPI(result);
	return ok ? 1 : 0;
}

LNLIB_EXPORT void LNLIB_NURBSSUR_non_uniform_scaling_query_size(const XYZW_C* control_points, int rows, int cols, double x_factor, double y_factor, double z_factor, const XYZ_C reference_point, XYZW_C** output_control_points, int* out_rows,int* out_cols)
{
	std::vector<std::vector<XYZW>> pts(rows, std::vector<XYZW>(cols));
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			pts[i][j] = ToXYZW(control_points[i * cols + j]);
		}
	}
	std::vector<std::vector<XYZW>> cps = NurbsSurface::NonUniformScaling(pts, x_factor, y_factor, z_factor, ToXYZ(reference_point));
	int r = cps.size();
	int c = cps[0].size();

	XYZW_C* result = new(std::nothrow) XYZW_C[rows * cols];
	if (!result) {
		return;
	}

	for (int i = 0; i < r; ++i) {
		for (int j = 0; j < c; ++j) {
			result[i * c + j] = FromXYZW(cps[i][j]);
		}
	}

	*output_control_points = result;
	*out_rows = r;
	*out_cols = c;
}

LNLIB_EXPORT void LNLIB_NURBSSUR_make_corner_fillet_surface(const LN_NurbsCurve_C arc, LN_NurbsSurface_C* out_surface)
{
	LN_NurbsSurface result;
	NurbsSurface::MakeCornerFilletSurface(FromCAPI(arc), result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_global_interpolation(const XYZ_C* points, int rows, int cols,
	int degree_u, int degree_v, LN_NurbsSurface_C* out_surface) 
{
	std::vector<std::vector<XYZ>> pts(rows, std::vector<XYZ>(cols));
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			pts[i][j] = ToXYZ(points[i * cols + j]);
		}
	}
	LN_NurbsSurface result;
	NurbsSurface::GlobalInterpolation(pts, degree_u, degree_v, result);
	*out_surface = ToCAPI(result);
}

int LNLIB_NURBSSUR_bicubic_local_interpolation(const XYZ_C* points, int rows, int cols, LN_NurbsSurface_C* out_surface) {
	std::vector<std::vector<XYZ>> pts(rows, std::vector<XYZ>(cols));
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			pts[i][j] = ToXYZ(points[i * cols + j]);
		}
	}
	LN_NurbsSurface result;
	bool ok = NurbsSurface::BicubicLocalInterpolation(pts, result);
	if (ok) *out_surface = ToCAPI(result);
	return ok ? 1 : 0;
}

int LNLIB_NURBSSUR_global_approximation(const XYZ_C* points, int rows, int cols,
	int degree_u, int degree_v, int ctrl_rows, int ctrl_cols, LN_NurbsSurface_C* out_surface) {
	std::vector<std::vector<XYZ>> pts(rows, std::vector<XYZ>(cols));
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			pts[i][j] = ToXYZ(points[i * cols + j]);
		}
	}
	LN_NurbsSurface result;
	bool ok = NurbsSurface::GlobalApproximation(pts, degree_u, degree_v, ctrl_rows, ctrl_cols, result);
	if (ok) *out_surface = ToCAPI(result);
	return ok ? 1 : 0;
}

int LNLIB_NURBSSUR_create_swung_surface(LN_NurbsCurve_C profile, LN_NurbsCurve_C trajectory,
	double scale, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	bool ok = NurbsSurface::CreateSwungSurface(FromCAPI(profile), FromCAPI(trajectory), scale, result);
	if (ok) *out_surface = ToCAPI(result);
	return ok ? 1 : 0;
}

void LNLIB_NURBSSUR_create_loft_surface(const LN_NurbsCurve_C* sections, int section_count,
	LN_NurbsSurface_C* out_surface, int custom_trajectory_degree,
	const double* custom_knots, int knot_count) {
	std::vector<LN_NurbsCurve> curves(section_count);
	for (int i = 0; i < section_count; ++i) {
		curves[i] = FromCAPI(sections[i]);
	}
	std::vector<double> knots;
	if (custom_knots && knot_count > 0) {
		knots.assign(custom_knots, custom_knots + knot_count);
	}
	LN_NurbsSurface result;
	NurbsSurface::CreateLoftSurface(curves, result, custom_trajectory_degree, knots);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_create_generalized_translational_sweep_surface(LN_NurbsCurve_C profile,
	LN_NurbsCurve_C trajectory, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::CreateGeneralizedTranslationalSweepSurface(FromCAPI(profile), FromCAPI(trajectory), result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_create_sweep_surface(LN_NurbsCurve_C profile, LN_NurbsCurve_C trajectory,
	int min_profiles, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::CreateSweepSurface(FromCAPI(profile), FromCAPI(trajectory), min_profiles, result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_create_sweep_surface_by_trajectory_degree(LN_NurbsCurve_C profile, LN_NurbsCurve_C trajectory,
	int min_profiles, int trajectory_degree, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::CreateSweepSurface(FromCAPI(profile), FromCAPI(trajectory), min_profiles, trajectory_degree, result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_create_gordon_surface(const LN_NurbsCurve_C* u_curves, int u_count,
	const LN_NurbsCurve_C* v_curves, int v_count,
	const XYZ_C* intersections, LN_NurbsSurface_C* out_surface) {
	std::vector<LN_NurbsCurve> ucs(u_count), vcs(v_count);
	for (int i = 0; i < u_count; ++i) ucs[i] = FromCAPI(u_curves[i]);
	for (int i = 0; i < v_count; ++i) vcs[i] = FromCAPI(v_curves[i]);
	std::vector<std::vector<XYZ>> ips(u_count, std::vector<XYZ>(v_count));
	for (int i = 0; i < u_count; ++i) {
		for (int j = 0; j < v_count; ++j) {
			ips[i][j] = ToXYZ(intersections[i * v_count + j]);
		}
	}
	LN_NurbsSurface result;
	NurbsSurface::CreateGordonSurface(ucs, vcs, ips, result);
	*out_surface = ToCAPI(result);
}

void LNLIB_NURBSSUR_create_coons_surface(LN_NurbsCurve_C left, LN_NurbsCurve_C bottom,
	LN_NurbsCurve_C right, LN_NurbsCurve_C top, LN_NurbsSurface_C* out_surface) {
	LN_NurbsSurface result;
	NurbsSurface::CreateCoonsSurface(FromCAPI(left), FromCAPI(bottom), FromCAPI(right), FromCAPI(top), result);
	*out_surface = ToCAPI(result);
}

double LNLIB_NURBSSUR_approximate_area(LN_NurbsSurface_C surface, int integrator_type) {
	return NurbsSurface::ApproximateArea(FromCAPI(surface), static_cast<IntegratorType>(integrator_type));
}

LNLIB_Mesh_C LNLIB_NURBSSUR_triangulate(LN_NurbsSurface_C surface, int resolution_u, int resolution_v, int use_delaunay) {
	LN_Mesh mesh = NurbsSurface::Triangulate(FromCAPI(surface), resolution_u, resolution_v, use_delaunay != 0);

	static thread_local std::vector<XYZ> g_vertices;
	static thread_local std::vector<int> g_faces_flat;
	static thread_local std::vector<UV> g_uvs;
	static thread_local std::vector<int> g_uv_indices;
	static thread_local std::vector<XYZ> g_normals;
	static thread_local std::vector<int> g_normal_indices;

	g_vertices = mesh.Vertices;

	g_faces_flat.clear();
	for (const auto& face : mesh.Faces) {
		for (int idx : face) {
			g_faces_flat.push_back(idx);
		}
	}

	g_uvs = mesh.UVs;

	if (!mesh.UVIndices.empty()) {
		g_uv_indices = mesh.UVIndices;
	}
	else {
		if (mesh.UVs.size() == mesh.Vertices.size()) {
			g_uv_indices = g_faces_flat;
		}
		else {
			g_uv_indices.clear();
		}
	}

	g_normals = mesh.Normals;

	if (!mesh.NormalIndices.empty()) {
		g_normal_indices = mesh.NormalIndices;
	}
	else {
		if (mesh.Normals.size() == mesh.Vertices.size()) {
			g_normal_indices = g_faces_flat;
		}
		else {
			g_normal_indices.clear();
		}
	}

	LNLIB_Mesh_C c{};
	c.vertices = reinterpret_cast<const XYZ_C*>(g_vertices.data());
	c.vertices_count = static_cast<int>(g_vertices.size());

	c.facets = g_faces_flat.data();
	c.facets_count = static_cast<int>(g_faces_flat.size());

	c.uvs = reinterpret_cast<const UV_C*>(g_uvs.data());
	c.uvs_count = static_cast<int>(g_uvs.size());

	c.uv_indices = g_uv_indices.data();
	c.uv_indices_count = static_cast<int>(g_uv_indices.size());

	c.normals = reinterpret_cast<const XYZ_C*>(g_normals.data());
	c.normals_count = static_cast<int>(g_normals.size());

	c.normal_indices = g_normal_indices.data();
	c.normal_indices_count = static_cast<int>(g_normal_indices.size());

	return c;
}

