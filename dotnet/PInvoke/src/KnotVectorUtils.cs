/*
 * Author:
 * 2025/11/30 - Yuqing Liang (BIMCoder Liang)
 * bim.frankliang@foxmail.com
 * 
 *
 * Use of this source code is governed by a LGPL-2.1 license that can be found in
 * the LICENSE file.
 */
using System;
using System.Runtime.InteropServices;
using static LNLibSharp.LNLibDefinitions;

namespace LNLibSharp
{
    public static partial class LNLibKnotVectorUtils
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_KV_get_continuity")]
        public static extern int GetContinuity(
        int degree,
        [In] double[] knot_vector,
        int knot_vector_count,
        double knot);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_KV_rescale")]
        public static extern void Rescale(
            [In] double[] knot_vector,
            int knot_vector_count,
            double min,
            double max,
            [Out] double[] result);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_KV_is_uniform")]
        public static extern int IsUniform(
            [In] double[] knot_vector,
            int knot_vector_count);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_KV_get_knot_multiplicity_map")]
        public static extern void GetKnotMultiplicityMap(
            [In] double[] knot_vector,
            int knot_vector_count,
            [In, Out] ref int out_size,
            [Out] double[] out_keys,
            [Out] int[] out_values);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_KV_get_internal_knot_multiplicity_map")]
        public static extern void GetInternalKnotMultiplicityMap(
            [In] double[] knot_vector,
            int knot_vector_count,
            [In, Out] ref int out_size,
            [Out] double[] out_keys,
            [Out] int[] out_values);
    }
}