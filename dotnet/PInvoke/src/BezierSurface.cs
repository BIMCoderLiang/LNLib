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
    public static partial class LNLibBezierSurface
    {
        [DllImport(
            LNLIB_CAPI_DLL, 
            EntryPoint = "LNLIB_BEZIERSURF_get_point_on_surface_by_deCasteljau",
            CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZ GetPointOnSurfaceByDeCasteljau(
            int degreeU,
            int degreeV,
            [In] XYZ[] control_points,
            int rows,
            int cols,
            UV uv);

        [DllImport(
           LNLIB_CAPI_DLL,
           EntryPoint = "LNLIB_BEZIERSURF_get_rational_point_on_surface_by_deCasteljau",
           CallingConvention = CallingConvention.Cdecl)]
        public static extern XYZW GetRationalPointOnSurfaceByDeCasteljau(
            int degreeU,
            int degreeV,
            [In] XYZW[] controlPoints,
            int rows,
            int cols,
            UV uv);
    }
}