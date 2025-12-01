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
    public static partial class LNLibAPI
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double polynomials_bernstein(int index, int degree, double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void polynomials_all_bernstein(int degree, double paramT, double[] out_array);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern double polynomials_horner_curve(
            int degree,
            IntPtr coefficients,
            int coeff_count,
            double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int polynomials_get_knot_span_index(
            int degree,
            IntPtr knot_vector,
            int knot_count,
            double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int polynomials_get_knot_multiplicity(
            IntPtr knot_vector,
            int knot_count,
            double paramT);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void polynomials_basis_functions(
            int span_index,
            int degree,
            IntPtr knot_vector,
            int knot_count,
            double paramT,
            double[] basis_functions);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern void polynomials_bezier_to_power_matrix(
            int degree,
            double[] out_matrix);
    }
}