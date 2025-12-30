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
    public static partial class LNLibBezierCurve
    {
        [DllImport(
            LNLIB_CAPI_DLL, 
            EntryPoint = "LNLIB_BEZIERCUR_get_point_on_curve_by_bernstein", 
            CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ GetPointOnCurveByBernstein(
            int degree,
            [In] XYZ[] controlPoints,
            int controlPointsCount,
            double paramT);

        [DllImport(
            LNLIB_CAPI_DLL,
            EntryPoint = "LNLIB_BEZIERCUR_get_rational_point_on_curve_by_bernstein",
            CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZW GetRationalPointOnCurveByBernstein(
            int degree,
            [In] XYZW[] controlPoints,
            int controlPointsCount,
            double paramT);

        [DllImport(
            LNLIB_CAPI_DLL,
            EntryPoint = "LNLIB_BEZIERCUR_get_point_on_curve_by_deCasteljau",
            CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ GetPointOnCurveByDeCasteljau(
            int degree,
            [In] XYZ[] controlPoints,
            int controlPointsCount,
            double paramT);

        [DllImport(
            LNLIB_CAPI_DLL,
            EntryPoint = "LNLIB_BEZIERCUR_get_rational_point_on_curve_by_deCasteljau",
            CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZW GetRationalPointOnCurveByDeCasteljau(
            int degree,
            [In] XYZW[] controlPoints,
            int controlPointsCount,
            double paramT);
    }
}
