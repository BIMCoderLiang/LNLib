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
        public static extern XYZ Create(double x, double y, double z);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_get_x")]
        public static extern double GetX(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_get_y")]
        public static extern double GetY(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_get_z")]
        public static extern double GetZ(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_set_x")]
        public static extern XYZ SetX(XYZ xyz, double x);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_set_y")]
        public static extern XYZ SetY(XYZ xyz, double y);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_set_z")]
        public static extern XYZ SetZ(XYZ xyz, double z);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_is_zero")]
        public static extern int IsZero(XYZ xyz, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_is_unit")]
        public static extern int IsUnit(XYZ xyz, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_is_almost_equal")]
        public static extern int IsAlmostEqual(XYZ a, XYZ b, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_length")]
        public static extern double Length(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_sqr_length")]
        public static extern double SqrLength(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_angle_to")]
        public static extern double AngleTo(XYZ a, XYZ b, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_normalize")]
        public static extern XYZ Normalize(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_add")]
        public static extern XYZ Add(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_subtract")]
        public static extern XYZ Subtract(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_negative")]
        public static extern XYZ Negative(XYZ xyz);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_dot_product")]
        public static extern double DotProduct(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZross_product")]
        public static extern XYZ CrossProduct(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_distance")]
        public static extern double Distance(XYZ a, XYZ b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_multiply")]
        public static extern XYZ Multiply(XYZ xyz, double scalar);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZ_divide")]
        public static extern XYZ Divide(XYZ xyz, double scalar);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_XYZreate_random_orthogonal")]
        public static extern XYZ CreateRandomOrthogonal(XYZ xyz);
    }
}