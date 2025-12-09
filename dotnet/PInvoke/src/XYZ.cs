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

    public static partial class LNLibXYZ
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZreate")]
        public static extern XYZ LNLIB_XYZreate(double x, double y, double z);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_get_x")]
        public static extern double LNLIB_XYZ_get_x(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_get_y")]
        public static extern double LNLIB_XYZ_get_y(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_get_z")]
        public static extern double LNLIB_XYZ_get_z(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_set_x")]
        public static extern XYZ LNLIB_XYZ_set_x(XYZ xyz, double x);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_set_y")]
        public static extern XYZ LNLIB_XYZ_set_y(XYZ xyz, double y);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_set_z")]
        public static extern XYZ LNLIB_XYZ_set_z(XYZ xyz, double z);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_is_zero")]
        public static extern int LNLIB_XYZ_is_zero(XYZ xyz, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_is_unit")]
        public static extern int LNLIB_XYZ_is_unit(XYZ xyz, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_is_almost_equal")]
        public static extern int LNLIB_XYZ_is_almost_equal(XYZ a, XYZ b, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_length")]
        public static extern double LNLIB_XYZ_length(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_sqr_length")]
        public static extern double LNLIB_XYZ_sqr_length(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_angle_to")]
        public static extern double LNLIB_XYZ_angle_to(XYZ a, XYZ b, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_normalize")]
        public static extern XYZ LNLIB_XYZ_normalize(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_add")]
        public static extern XYZ LNLIB_XYZ_add(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_subtract")]
        public static extern XYZ LNLIB_XYZ_subtract(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_negative")]
        public static extern XYZ LNLIB_XYZ_negative(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_dot_product")]
        public static extern double LNLIB_XYZ_dot_product(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZross_product")]
        public static extern XYZ LNLIB_XYZross_product(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_distance")]
        public static extern double LNLIB_XYZ_distance(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_multiply")]
        public static extern XYZ LNLIB_XYZ_multiply(XYZ xyz, double scalar);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_divide")]
        public static extern XYZ LNLIB_XYZ_divide(XYZ xyz, double scalar);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZreate_random_orthogonal")]
        public static extern XYZ LNLIB_XYZreate_random_orthogonal(XYZ xyz);
    }
}