/*
 * Author:
 * 2025/12/09 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

using System.Runtime.InteropServices;
using static LNLibSharp.LNLibDefinitions;

namespace LNLibSharp
{
    public static partial class LNLibInterpolation
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_get_total_chord_length")]
        public static extern double GetTotalChordLength(
        [In] XYZ[] through_points,
        int through_points_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_get_chord_parameterization")]
        public static extern void GetChordParameterization(
            [In] XYZ[] through_points,
            int through_points_count,
            [Out] double[] out_params,
            [In, Out] ref int out_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_get_centripetal_length")]
        public static extern double GetCentripetalLength(
            [In] XYZ[] through_points,
            int through_points_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_get_centripetal_parameterization")]
        public static extern void GetCentripetalParameterization(
            [In] XYZ[] through_points,
            int through_points_count,
            [Out] double[] out_params,
            [In, Out] ref int out_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_average_knot_vector")]
        public static extern void AverageKnotVector(
            int degree,
            [In] double[] parameters,
            int params_count,
            [Out] double[] out_knot_vector,
            [In, Out] ref int out_knot_vector_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_get_surface_mesh_parameterization")]
        public static extern int GetSurfaceMeshParameterization(
            [In] XYZ[] through_points,
            int rows,
            int cols,
            [Out] double[] params_u,
            [In, Out] ref int out_params_u_size,
            [Out] double[] params_v,
            [In, Out] ref int out_params_v_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_compute_tangent_by_three_points")]
        public static extern void ComputeTangentByThreePoints(
            [In] XYZ[] through_points,
            int through_points_count,
            [Out] XYZ[] out_tangents,
            [In, Out] ref int out_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_compute_tangent_by_five_points")]
        public static extern int ComputeTangentByFivePoints(
            [In] XYZ[] through_points,
            int through_points_count,
            [Out] XYZ[] out_tangents,
            [In, Out] ref int out_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_computer_weight_for_rational_quadratic_interpolation")]
        public static extern int ComputerWeightForRationalQuadraticInterpolation(
            XYZ start_point,
            XYZ middle_control_point,
            XYZ end_point,
            out double out_weight);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_INTERPOLATION_compute_knot_vector")]
        public static extern void ComputeKnotVector(
            int degree,
            int control_points_count,
            [In] double[] parameters,
            int params_count,
            [Out] double[] out_knot_vector,
            [In, Out] ref int out_knot_vector_size);
    }
}