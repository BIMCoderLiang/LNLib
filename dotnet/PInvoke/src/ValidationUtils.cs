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
using static LNLibSharp.LNLibDefinitions;

namespace LNLibSharp
{
    public static partial class LNLibValidationUtils
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_VALID_is_valid_knotVector")]
        public static extern int IsValidKnotVector(
        [In] double[] knot_vector,
        int knot_vector_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_VALID_is_valid_bezier")]
        public static extern int IsValidBezier(
            int degree,
            int control_points_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_VALID_is_valid_bspline")]
        public static extern int IsValidBspline(
            int degree,
            int knot_vector_count,
            int control_points_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_VALID_is_valid_nurbs")]
        public static extern int IsValidNurbs(
            int degree,
            int knot_vector_count,
            int control_points_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_VALID_compute_curve_modify_tolerance")]
        public static extern double ComputeCurveModifyTolerance(
            [In] XYZW[] control_points,
            int control_points_count);
    }
}