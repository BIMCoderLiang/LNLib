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
    public static extern XYZ bezier_surface_get_point_by_de_casteljau(
        int degree_u,
        int degree_v,
        IntPtr control_points,
        int num_u,
        int num_v,
        UV uv);
}
}