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
    public static partial class LNLibPolynomials
    {
        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_horner")]
        public static extern double LNLIB_POLY_horner(
        int degree,
        [In] double[] coefficients,
        int coeff_count,
        double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_bernstein")]
        public static extern double LNLIB_POLY_bernstein(
            int index,
            int degree,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_all_bernstein")]
        public static extern void LNLIB_POLY_all_bernstein(
            int degree,
            double param_t,
            [Out] double[] out_array);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_horner_uv")]
        public static extern void LNLIB_POLY_horner_uv(
            int degree_u,
            int degree_v,
            [In] double[] coefficients,
            UV uv,
            [Out] double[] out_horner);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_get_knot_multiplicity")]
        public static extern int LNLIB_POLY_get_knot_multiplicity(
            [In] double[] knot_vector,
            int knot_vector_count,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_get_knot_span_index")]
        public static extern int LNLIB_POLY_get_knot_span_index(
            int degree,
            [In] double[] knot_vector,
            int knot_vector_count,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_basis_functions")]
        public static extern void LNLIB_POLY_basis_functions(
            int span_index,
            int degree,
            [In] double[] knot_vector,
            int knot_vector_count,
            double param_t,
            [Out] double[] out_basis_functions);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_basis_functions_derivatives")]
        public static extern void LNLIB_POLY_basis_functions_derivatives(
            int span_index,
            int degree,
            int derivative,
            [In] double[] knot_vector,
            int knot_vector_count,
            double param_t,
            out IntPtr out_derivatives,
            out int out_rows,
            out int out_cols);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_basis_functions_first_derivative")]
        public static extern void LNLIB_POLY_basis_functions_first_derivative(
            int span_index,
            int degree,
            [In] double[] knot_vector,
            int knot_vector_count,
            double param_t,
            [Out] double[] out_derivatives);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_one_basis_function")]
        public static extern double LNLIB_POLY_one_basis_function(
            int span_index,
            int degree,
            [In] double[] knot_vector,
            int knot_vector_count,
            double param_t);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_one_basis_function_derivative")]
        public static extern void LNLIB_POLY_one_basis_function_derivative(
            int span_index,
            int degree,
            int derivative,
            [In] double[] knot_vector,
            int knot_vector_count,
            double param_t,
            out IntPtr out_derivatives,
            out int out_size);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_all_basis_functions")]
        public static extern void LNLIB_POLY_all_basis_functions(
            int span_index,
            int degree,
            [In] double[] knot_vector,
            int knot_vector_count,
            double knot,
            out IntPtr out_basis_functions,
            out int out_rows,
            out int out_cols);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_bezier_to_power_matrix")]
        public static extern void LNLIB_POLY_bezier_to_power_matrix(
            int degree,
            out IntPtr out_matrix,
            out int out_rows,
            out int out_cols);

        [DllImport(LNLIB_CAPI_DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "LNLIB_POLY_power_to_bezier_matrix")]
        public static extern void LNLIB_POLY_power_to_bezier_matrix(
            int degree,
            [In] double[] matrix,
            out IntPtr out_matrix,
            out int out_rows,
            out int out_cols);
    }
}