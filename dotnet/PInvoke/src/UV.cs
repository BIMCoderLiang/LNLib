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
using static LNLibSharp.LNLibDefinitions;

namespace LNLibSharp
{
    [StructLayout(LayoutKind.Sequential)]
    public struct UV
    {
        public double u, v;
    }

    public static partial class LNLibUV
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UVreate")]
        public static extern UV Create(double u, double v);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_get_u")]
        public static extern double GetU(UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_get_v")]
        public static extern double GetV(UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_set_u")]
        public static extern UV SetU(UV uv, double u);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_set_v")]
        public static extern UV SetV(UV uv, double v);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_is_zero")]
        public static extern int IsZero(UV uv, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_is_unit")]
        public static extern int IsUnit(UV uv, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_is_almost_equal")]
        public static extern int IsAlmostEqual(UV a, UV b, double epsilon);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_length")]
        public static extern double Length(UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_sqr_length")]
        public static extern double SqrLength(UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_normalize")]
        public static extern UV Normalize(UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_add")]
        public static extern UV Add(UV a, UV b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_subtract")]
        public static extern UV Subtract(UV a, UV b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_negative")]
        public static extern UV Negative(UV uv);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_dot_product")]
        public static extern double DotProduct(UV a, UV b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UVross_product")]
        public static extern double CrossProduct(UV a, UV b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_distance")]
        public static extern double Distance(UV a, UV b);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_scale")]
        public static extern UV Scale(UV uv, double factor);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_UV_divide")]
        public static extern UV Divide(UV uv, double divisor);
    }
}