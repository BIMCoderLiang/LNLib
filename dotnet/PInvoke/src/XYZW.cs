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
    [StructLayout(LayoutKind.Sequential)]
    public struct XYZW
    {
        public double wx, wy, wz, w;
    }

    public static partial class LNLibXYZW
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_create")]
        public static extern XYZW Create(double wx, double wy, double wz, double w);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_create_from_xyz")]
        public static extern XYZW CreateFromXYZ(XYZ xyz, double w);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_get_wx")]
        public static extern double GetWX(XYZW xyzw);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_get_wy")]
        public static extern double GetWY(XYZW xyzw);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_get_wz")]
        public static extern double GetWZ(XYZW xyzw);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_get_w")]
        public static extern double GetW(XYZW xyzw);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_set_w")]
        public static extern XYZW SetW(XYZW xyzw, double w);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_to_xyz")]
        public static extern XYZ ToXYZ(XYZW xyzw, int divideWeight);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_is_almost_equal")]
        public static extern int IsAlmostEqual(XYZW a, XYZW b, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_distance")]
        public static extern double Distance(XYZW a, XYZW b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_add")]
        public static extern XYZW Add(XYZW a, XYZW b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_subtract")]
        public static extern XYZW Subtract(XYZW a, XYZW b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_multiply")]
        public static extern XYZW Multiply(XYZW a, double scalar);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZW_divide")]
        public static extern XYZW Divide(XYZW a, double scalar);
    }
}