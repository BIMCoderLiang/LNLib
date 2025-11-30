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
    public static extern XYZ bezier_curve_get_point_by_bernstein(
        int degree,
        IntPtr control_points,
        int control_points_count,
        double paramT);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZ bezier_curve_get_point_by_de_casteljau(
        int degree,
        IntPtr control_points,
        int control_points_count,
        double paramT);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZW bezier_curve_get_point_by_bernstein_rational(
        int degree,
        IntPtr control_points,
        int control_points_count,
        double paramT);
}
} 
