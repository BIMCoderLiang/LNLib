/*
 * Author:
 * 2025/11/30 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */

namespace LNLibSharp
{
public static partial class LNLibAPI
{
    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int validation_utils_is_valid_knot_vector(
        IntPtr knot_vector,
        int count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int validation_utils_is_valid_bezier(
        int degree,
        int control_points_count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int validation_utils_is_valid_bspline(
        int degree,
        int knot_count,
        int cp_count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int validation_utils_is_valid_nurbs(
        int degree,
        int knot_count,
        int weighted_cp_count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double validation_utils_compute_curve_modify_tolerance(
        IntPtr control_points,
        int count);
}
}