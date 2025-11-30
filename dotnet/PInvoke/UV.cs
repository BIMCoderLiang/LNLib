/*
 * Author:
 * 2025/11/25 - Yuqing Liang (BIMCoder Liang)
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
public struct UV
{
    public double u, v;
}

public static partial class LNLibAPI
{
    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV uv_create(double u, double v);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double uv_get_u(UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double uv_get_v(UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV uv_add(UV a, UV b);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV uv_subtract(UV a, UV b);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV uv_negative(UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV uv_normalize(UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV uv_scale(UV uv, double factor);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern UV uv_divide(UV uv, double divisor);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double uv_length(UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double uv_sqr_length(UV uv);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double uv_distance(UV a, UV b);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int uv_is_zero(UV uv, double epsilon);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int uv_is_unit(UV uv, double epsilon);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int uv_is_almost_equal(UV a, UV b, double epsilon);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double uv_dot(UV a, UV b);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double uv_cross(UV a, UV b);
}
}