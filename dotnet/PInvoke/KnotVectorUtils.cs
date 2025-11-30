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
    public static extern int knot_vector_utils_get_continuity(
        int degree,
        IntPtr knot_vector,
        int knot_vector_count,
        double knot);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void knot_vector_utils_rescale(
        IntPtr knot_vector,
        int knot_vector_count,
        double min_val,
        double max_val,
        double[] out_rescaled_knot_vector);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int knot_vector_utils_is_uniform(
        IntPtr knot_vector,
        int knot_vector_count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int knot_vector_utils_get_knot_multiplicity_map_size(
        IntPtr knot_vector,
        int knot_vector_count);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void knot_vector_utils_get_knot_multiplicity_map(
        IntPtr knot_vector,
        int knot_vector_count,
        double[] out_unique_knots,
        int[] out_multiplicities);
}
}