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
public struct XYZW
{
    public double wx, wy, wz, w;
}

public static partial class LNLibAPI
{
    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZW xyzw_create(double wx, double wy, double wz, double w);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZW xyzw_create_from_xyz(XYZ xyz, double w);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZ xyzw_to_xyz(XYZW v, int divideWeight);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZW xyzw_add(XYZW a, XYZW b);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZW xyzw_multiply(XYZW a, double scalar);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern XYZW xyzw_divide(XYZW a, double scalar);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double xyzw_distance(XYZW a, XYZW b);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double xyzw_get_wx(XYZW v);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double xyzw_get_wy(XYZW v);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double xyzw_get_wz(XYZW v);

    [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern double xyzw_get_w(XYZW v);
}
}