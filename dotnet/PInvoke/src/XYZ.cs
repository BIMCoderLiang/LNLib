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
    public struct XYZ
    {
        public double x, y, z;
    }

    public static partial class LNLibAPI
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_create(double x, double y, double z);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_zero();

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_add(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_subtract(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_negative(XYZ a);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_multiply(XYZ a, double scalar);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_divide(XYZ a, double scalar);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double xyz_length(XYZ v);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double xyz_sqr_length(XYZ v);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int xyz_is_zero(XYZ v, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int xyz_is_unit(XYZ v, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_normalize(XYZ v);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double xyz_dot(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ xyz_cross(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double xyz_distance(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int xyz_equals(XYZ a, XYZ b);
    }
}